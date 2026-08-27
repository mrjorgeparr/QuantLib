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

`Greeks` layers delta, gamma, vega, rho and theta on top of any of them by finite differences.

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
    //                  S0      K      r    sigma   T   isPut
    EuropeanOption euro(100.0, 100.0, 0.05, 0.20, 1.0, true);
    AmericanOption amer(100.0, 100.0, 0.05, 0.20, 1.0, true);

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
│   ├── european_option.hpp  EuropeanOption : BaseOption  (exercise-style tag type)
│   └── american_option.hpp  AmericanOption : BaseOption  (exercise-style tag type)
├── greeks/
│   └── greeks.hpp          Greeks — finite-difference sensitivities
├── math/
│   └── lsmc_utils.hpp       MathUtils — Laguerre basis, least squares, Gaussian elimination
└── pricers/
    └── base_pricer.hpp      BasePricer — the pricing engine (paths + all three methods)

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

## Greeks

`Greeks` computes delta, gamma, vega, rho and theta by finite differences, driving whatever pricer you hand it.

```cpp
#include "greeks/greeks.hpp"

EuropeanOption euro(100.0, 100.0, 0.05, 0.20, 1.0, true);
BasePricer pricer(50, 200000);

Greeks greeks(euro, pricer);        // works with any instrument the pricer prices
greeks.setSeed(20260827u);          // optional: makes the run reproducible

GreekSet g = greeks.all();          // or greeks.delta(), greeks.vega(), ...
printf("delta %.4f  gamma %.5f  vega %.2f\n", g.delta, g.gamma, g.vega);
```

Validated against closed-form Black-Scholes (S=K=100, r=5%, sigma=20%, T=1, 200k paths):

| | Greeks | Black-Scholes | error |
| --- | ---: | ---: | ---: |
| price | 5.60079 | 5.57353 | +0.49% |
| delta | -0.36480 | -0.36317 | -0.45% |
| gamma | 0.01874 | 0.01876 | -0.11% |
| vega | 37.68755 | 37.52403 | +0.44% |
| rho | -42.07588 | -41.89046 | -0.44% |
| theta | -1.66527 | -1.65788 | -0.45% |

### Only the constructor is a template

`BasePricer::price` is a set of compile-time **overloads** — it dispatches on the *static* type of its argument. A member of type `BaseOption` would both slice the object and match no overload.

`Greeks` solves that without becoming a template itself. The contract lives in a plain `BaseOption option_` member, and bumped scenarios are built with `BaseOption` and `Ticker` — `setTicker`, `setK`, `setR`, `setT`, and their validation. Nothing about the contract is re-modelled.

Only the *dispatch* needs the concrete type, so only the constructor is templated: it seals a copy of the instrument in a `std::function<double(const BaseOption&)>` that transplants bumped terms onto a correctly-typed copy and calls `price()`. That one function is the entire cost of the overload set; every finite-difference method is ordinary code in `src/greeks/greeks.cpp`.

Two payoffs. `Greeks` is one plain type, so `std::vector<Greeks>` holds European and American positions side by side. And because the lambda copies the *concrete* option each call, an instrument that later carries its own state (a barrier level, an averaging schedule) keeps it — a sliced `BaseOption` would not. The cost is one indirect call per pricing, nothing against a 200k-path Monte Carlo.

### Common random numbers are not optional

Each Greek is a difference of two Monte Carlo estimates divided by a small bump. If the two legs use independent draws, their sampling error does not cancel and the division amplifies what remains. Measured on this library, delta with an *absolute* 0.01 bump and independent draws:

| | mean | std dev | observed range |
| --- | ---: | ---: | --- |
| Independent draws | -0.207 | 2.027 | [-4.09, +2.29] |
| Common random numbers | -0.3631 | 0.0014 | [-0.365, -0.361] |

A put delta is confined to [-1, 0]; the first row gets the *sign* wrong about a third of the time. `Greeks` therefore pins the seed on both legs by default (`setCommonRandomNumbers(false)` disables it, which is useful only for demonstrating the above). This is what `BasePricer::setSeed` exists for. A `SeedGuard` restores the pricer seed state on scope exit, so computing Greeks never leaves a caller pricer reconfigured.

### Bump sizes and units

Spot bumps are **relative** (a fraction of S0); the rest are absolute. Gamma gets a larger default bump than delta because the second difference divides by `h^2`, which amplifies noise quadratically — at a 1%-of-spot bump gamma carries ~16% noise even with CRN, versus ~0.7% at 5%.

| Greek | default bump | units of the result |
| --- | --- | --- |
| delta | 1% of spot | per unit of spot |
| gamma | 5% of spot | per unit of spot squared |
| vega | 0.01 vol | per 1.0 of vol (divide by 100 for per vol point) |
| rho | 0.001 | per 1.0 of rate (divide by 10000 for per bp) |
| theta | 1/365 year | per year (divide by 365 for per day) |

Override with `GreekBumps` via the three-argument constructor or `setBumps`. Vega and theta fall back to a one-sided difference when the bump would push sigma or the maturity through zero.

`all()` costs 10 pricer calls; the individual accessors cost 1-3 each. With the LSMC hot spot in place, a full American `all()` at 100k paths is several minutes — worth fixing gap 2 first if you plan to work with American Greeks.

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

1. **`AmericanOption` carries an unused RNG.** `gen` and `norm_dist` are leftovers from when the instrument priced itself; path generation now lives in `BasePricer::generateGBMPaths`. They are marked `TODO(unused)` in the header. Removing them also lets `<random>` drop out of the header.
2. **LSMC has a quadratic hot spot.** The "roll forward the out-of-the-money paths" step does a linear scan of `inMoneyIndices` for every path (`src/pricers/base_pricer.cpp`, the `for (int index : inMoneyIndices)` loop) — O(paths x in-the-money paths) per time step. Replacing it with a `std::vector<char>` membership mask makes it O(paths) and is a large speedup at high path counts.
3. **`BasePricer::price` overloads are not virtual**, despite the `virtual ~BasePricer()`. The polymorphism the destructor implies is not actually wired up. Decide: plain overloads on a concrete class, or a virtual `price(const BaseOption&)` interface.
4. **`/openmp` is set twice on MSVC** — once in `CMAKE_CXX_FLAGS_RELEASE`, once via `OpenMP::OpenMP_CXX`. Harmless, but drop it from the manual flags.
5. **No tests, no benchmarks, no CI.** `Greeks` can be validated against closed-form Black-Scholes (it currently agrees to within 0.5% on every Greek); that is the obvious first test to commit.
6. **Gamma on short-dated or far out-of-the-money options** may need a bump size other than the 5%-of-spot default — the bias/variance tradeoff moves with the contract.
7. **Everything is in the global namespace.** `Ticker` in particular is a name that will collide. Worth fixing at the same time as the rename.

### Recently resolved

- Orphaned `EuropeanPricer` / `LSMC` declarations moved out of `include/pricers/` into `deprecated/`, so including a header with no definition is no longer possible.
- `build/` untracked from git and a `.gitignore` added.
- `N`, `M` and `dt` removed from `EuropeanOption` / `AmericanOption`. They were stored but ignored by `BasePricer`, so those constructor arguments silently did nothing. Discretisation now lives only on the pricer.
- `BasePricer::setSeed` / `clearSeed` added, so a run can be made reproducible. This was a prerequisite for the Greeks, not just a convenience — see below.
- Greeks implemented in `Greeks` (`include/greeks/greeks.hpp`, `src/greeks/greeks.cpp`).

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
