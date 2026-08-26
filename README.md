# MyQuantLib

> A small, dependency-free C++11 static library for Monte Carlo option pricing, parallelised with OpenMP.

**Status:** early / working prototype. The name `MyQuantLib` is a placeholder and will change — see [Renaming the library](#renaming-the-library).

---

## What it does

One engine, `BasePricer`, generates geometric Brownian motion paths and prices three payoff styles off them:

| Method | Instrument | Technique |
| --- | --- | --- |
| `price(const EuropeanOption&)` | European call/put | Terminal-payoff Monte Carlo |
| `price(const AmericanOption&)` | American call/put | Longstaff–Schwartz (LSMC) least-squares regression |
| `priceArithmeticAsian(const BaseOption&)` | Arithmetic-average Asian | Path-average Monte Carlo |

No external dependencies. No BLAS, no Boost — the least-squares solve is hand-rolled (normal equations + Gaussian elimination with partial pivoting) over four Laguerre basis functions.

---

## Quick start

### Build (CMake — recommended)

```bash
cmake -S . -B build
cmake --build build --config Release
```

Produces `build/Release/MyQuantLib.lib` (MSVC) or `build/libMyQuantLib.a` (GCC/Clang).

### Build (Makefile — Unix / MSYS2 / Git Bash)

```bash
make
```

Produces `libMyQuantLib.a` in the project root.

### Use it

```cpp
#include "instruments/european_option.hpp"
#include "instruments/american_option.hpp"
#include "pricers/base_pricer.hpp"
#include <cstdio>

int main() {
    //                  S0      K      r    sigma   T   N     M     isPut
    EuropeanOption euro(100.0, 100.0, 0.05, 0.20, 1.0, 50, 100000, true);
    AmericanOption amer(100.0, 100.0, 0.05, 0.20, 1.0, 50, 100000, true);

    BasePricer pricer(/*timeSteps=*/50, /*paths=*/100000);

    printf("European put: %.4f\n", pricer.price(euro));                 // ~5.55  (Black-Scholes: 5.5735)
    printf("American put: %.4f\n", pricer.price(amer));                 // ~6.28  (early-exercise premium)
    printf("Asian put   : %.4f\n", pricer.priceArithmeticAsian(euro));  // ~3.38
    return 0;
}
```

Compile against the library:

```bash
g++ -std=c++11 -O3 -fopenmp -Iinclude demo.cpp -L. -lMyQuantLib -o demo
```

---

## Architecture

```
include/
├── core/
│   ├── constants.hpp        OptionConstants: default steps/paths, FD bump sizes, EPSILON
│   └── ticker.hpp           Ticker — the underlying asset: (S0, sigma)
├── instruments/
│   ├── base_option.hpp      BaseOption — the contract: Ticker + K, r, T, isPut, payoff()
│   ├── european_option.hpp  EuropeanOption : BaseOption
│   └── american_option.hpp  AmericanOption : BaseOption
├── math/
│   └── lsmc_utils.hpp       MathUtils — Laguerre basis, least squares, Gaussian elimination
└── pricers/
    ├── base_pricer.hpp      BasePricer — the pricing engine (paths + all three methods)
    ├── european_pricer.hpp  (declaration only — see Known gaps)
    └── lsmc.hpp             (declaration only — see Known gaps)

src/                         mirrors include/, one .cpp per header that needs one
deprecated/pricers/          pre-refactor pricers, NOT compiled
```

The design separates three concerns deliberately:

1. **`Ticker`** — what the *asset* is (spot, volatility). Reusable across contracts.
2. **`BaseOption`** — what the *contract* is (strike, rate, maturity, put/call) and how it pays off. `payoff()` is the only place call/put asymmetry lives, so every pricing method supports both for free.
3. **`BasePricer`** — *how* you value it (how many paths, how many steps). Discretisation is a property of the numerical method, not of the contract, so it lives here.

### Path generation

`generateGBMPaths` uses the exact GBM solution rather than an Euler step:

```
S(t+dt) = S(t) * exp[ (r - 0.5*sigma^2)*dt + sigma*sqrt(dt)*Z ],   Z ~ N(0,1)
```

so there is **no time-discretisation bias for European payoffs** — only Monte Carlo sampling error. Path-dependent payoffs (Asian, American) still carry discretisation error from monitoring at `N` discrete dates.

Each OpenMP thread gets its own `std::mt19937`, seeded from a shared draw plus `omp_get_thread_num()`. This keeps generation thread-safe, but it also means results are **not reproducible across runs or thread counts**.

### Longstaff–Schwartz

Standard backward induction:

1. Set terminal cash flows to the terminal payoff.
2. Walk backwards from `N-1` to `1`. At each date, take only the in-the-money paths.
3. Regress discounted future cash flow on `{1, L0(x), L1(x), L2(x)}` where `x = S/K` (moneyness — scaling keeps the exponential basis numerically well behaved).
4. Exercise where immediate payoff beats the fitted continuation value; otherwise roll the cash flow forward.
5. Discount the surviving cash flows at `t=1` back to today.

The regression is skipped when fewer than `BASIS_FUNCTIONS_COUNT` paths are in the money — with 4 unknowns and 4 or fewer observations the normal-equations matrix is singular.

---

## Build system

Two independent build paths that produce the same artifact. **CMake is the primary one**; the Makefile is a GNU-make fallback.

### CMakeLists.txt

| Directive | Why it's there |
| --- | --- |
| `cmake_minimum_required(VERSION 3.10)` | Floor for the `OpenMP::OpenMP_CXX` imported target. |
| `set(CMAKE_CXX_STANDARD 11)` + `..._REQUIRED ON` | C++11, and *fail* rather than silently downgrade. |
| `find_package(OpenMP REQUIRED)` | Locates the compiler's OpenMP flags/libs; configure fails without it. |
| `if(MSVC) ... else()` | MSVC and GCC/Clang have completely different flag spellings (`/O2` vs `-O3`). |
| `add_compile_options(/wd4267)` | Silences MSVC's `size_t` → `int` narrowing warning (from the `.size()` casts). |
| `file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")` | Auto-discovers sources. `CONFIGURE_DEPENDS` makes CMake re-check the glob on every build, so new files are picked up without manually re-running configure. Note it globs **only `src/`** — which is why `deprecated/` is excluded for free. |
| `add_library(... STATIC ...)` | A static archive: not an executable, not a DLL. |
| `target_link_libraries(... PUBLIC OpenMP::OpenMP_CXX)` | `PUBLIC` propagates OpenMP's flags and runtime to anything that links this library — consumers don't have to remember `-fopenmp`. |
| `target_include_directories(... PUBLIC include)` | Same idea for headers: consumers get `include/` on their search path automatically. |

**Multi-config caveat.** The Visual Studio generator is *multi-config*: `CMAKE_BUILD_TYPE` at configure time is ignored, and you choose the configuration at build time with `--config Release`. Single-config generators (Makefiles, Ninja) are the reverse — `-DCMAKE_BUILD_TYPE=Release` at configure time, no `--config`. Getting this backwards is the usual cause of "why did I get an unoptimised build".

### Makefile

```make
SOURCES = $(shell find $(SRCDIR) -name '*.cpp')
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
```

`find` collects the sources; `patsubst` rewrites `src/pricers/base_pricer.cpp` into `obj/pricers/base_pricer.o`, mirroring the source tree under `obj/`.

```make
$(TARGET): $(OBJECTS)
	ar rcs $@ $^
```

The library depends on every object file. `ar rcs` creates the archive (`r` insert, `c` create quietly, `s` write the symbol index). `$@` is the target, `$^` all prerequisites.

```make
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
```

A pattern rule: any `.o` is built from the matching `.cpp`. `mkdir -p $(dir $@)` creates the nested output directory, since `patsubst` invents subdirectories that don't exist yet. `$<` is the first prerequisite (the `.cpp`), `$@` the output. Make recompiles an object only when its `.cpp` is newer — that's the incremental build.

```make
debug: CXXFLAGS += -g -DDEBUG
debug: $(TARGET)
```

A *target-specific variable*: `-g -DDEBUG` is appended to `CXXFLAGS` only for the `debug` target and everything it pulls in.

> **`make debug` straight after `make` does nothing.** The objects built in Release are already newer than their sources, so Make considers them up to date and never recompiles with the debug flags. Run `make clean debug`.

> **The Makefile needs a POSIX shell.** `find` and `mkdir -p` mean Git Bash / MSYS2 / WSL / Linux / macOS. Under `cmd.exe`, Windows' own `find.exe` is a text-search tool and the build silently produces nothing. Use CMake on native Windows.

### Two differences worth knowing

- **Header dependencies.** CMake tracks them (edit `base_option.hpp`, everything that includes it rebuilds). The Makefile does **not** — it only knows that a `.o` depends on its `.cpp`. After editing a header, run `make clean`.
- **Warnings.** The Makefile compiles with `-Wall -Wextra`; the CMake build does not. The Makefile is the stricter of the two.

---

## Known gaps

An honest list of what is unfinished, roughly in priority order.

1. **`EuropeanPricer` and `LSMC` are declared but not defined.** `include/pricers/european_pricer.hpp` and `include/pricers/lsmc.hpp` still ship, but their implementations were moved to `deprecated/pricers/` during the refactor that folded everything into `BasePricer` — and `deprecated/` is not compiled. Including either header and instantiating the class gives an undefined-reference link error. (`deprecated/pricers/european_pricer.cpp` is empty.) Either delete both headers or reinstate them as thin wrappers.
2. **`build/` is committed to git** — roughly 100 MSVC/MSBuild artifacts, including stale objects from before the refactor. Add a `.gitignore` (`build/`, `obj/`, `*.a`, `*.lib`, `*.obj`, `*.o`) and `git rm -r --cached build`.
3. **Discretisation parameters are duplicated and inert.** `EuropeanOption` and `AmericanOption` each store `N`, `M`, `dt`, but `BasePricer` ignores them and uses its own `timeSteps_`/`paths_`. Those constructor arguments silently do nothing — a real trap for a caller. Drop them from the instruments (they belong to the pricer) or have the pricer read them.
4. **`AmericanOption` carries an unused RNG.** `gen` and `norm_dist` are leftovers from when the instrument priced itself; path generation now lives in `BasePricer`.
5. **LSMC has a quadratic hot spot.** The "roll forward the out-of-the-money paths" step does a linear scan of `inMoneyIndices` for every path (`src/pricers/base_pricer.cpp`, the `for (int index : inMoneyIndices)` loop) — O(paths × in-the-money paths) per time step. Replacing it with a `std::vector<char>` membership mask is a large speedup at high path counts.
6. **`BasePricer::price` overloads are not virtual**, despite the `virtual ~BasePricer()`. The polymorphism the destructor implies is not actually wired up. Decide: plain overloads on a concrete class, or a virtual `price(const BaseOption&)` interface.
7. **`/openmp` is set twice on MSVC** — once in `CMAKE_CXX_FLAGS_RELEASE`, once via `OpenMP::OpenMP_CXX`. Harmless, but drop it from the manual flags.
8. **No Greeks.** `constants.hpp` defines finite-difference bump sizes (`FINITE_DIFF_BUMP_PRICE`, `..._VOL`, `..._RATE`, `..._TIME`) but nothing uses them — delta/vega/rho/theta were planned and never written.
9. **No tests, no benchmarks, no CI.** The European price can be validated against closed-form Black–Scholes; that is the obvious first test.
10. **Results are not reproducible.** There is no way to fix the seed. A seeded `BasePricer` constructor would make runs deterministic and testable.

---

## Renaming the library

The name appears in exactly **five** places, all in the two build files. Nothing in `include/` or `src/` hard-codes it (the include guards are generic), so the rename is cheap:

| File | Line | Current |
| --- | --- | --- |
| `CMakeLists.txt` | 2 | `project(MyQuantLib)` |
| `CMakeLists.txt` | 29 | `add_library(MyQuantLib STATIC ${SOURCES})` |
| `CMakeLists.txt` | 33 | `target_link_libraries(MyQuantLib PUBLIC ...)` |
| `CMakeLists.txt` | 36 | `target_include_directories(MyQuantLib PUBLIC include)` |
| `Makefile` | 5 | `TARGET = libMyQuantLib.a` |

Delete `build/` before renaming — a stale `CMakeCache.txt` pins the old project name and source directory.

If the code ever gets a namespace (worth doing before this is used anywhere real), that is the moment to add it: everything currently sits in the global namespace, and `Ticker` in particular is a name that will collide.

---

## Requirements

- A C++11 compiler with OpenMP: MSVC 2019+, GCC 7+, or Clang with `libomp`
- CMake 3.10+ *(or)* GNU Make with a POSIX shell
- No third-party libraries
