# Ito

> A small, dependency-free C++11 library for Monte Carlo option pricing and Greeks, parallelised with OpenMP.

Named for the Itô processes it simulates — geometric Brownian motion is the canonical one.

**Status:** working prototype. Prices and Greeks are validated against closed-form Black–Scholes; there is no test suite committed yet.

---

## What it does

One engine, `BasePricer`, generates GBM paths and values three payoff styles from them:

| Call | Instrument | Method |
| --- | --- | --- |
| `price(const EuropeanOption&)` | European call/put | Terminal-payoff Monte Carlo |
| `price(const AmericanOption&)` | American call/put | Longstaff–Schwartz least-squares regression |
| `priceArithmeticAsian(const BaseOption&)` | Arithmetic-average Asian | Path-average Monte Carlo |

`Greeks` layers delta, gamma, vega, rho and theta on top of any of them by finite differences.

No third-party libraries. The least-squares solve is hand-rolled — normal equations plus Gaussian elimination with partial pivoting, over four Laguerre basis functions.

---

## Quick start

### Build

```bash
cmake -S . -B build
cmake --build build --config Release
```

Produces `build/Release/Ito.lib` (MSVC) or `build/libIto.a` (GCC/Clang).

A GNU Makefile is also provided for Unix, MSYS2 or Git Bash:

```bash
make
```

### Price something

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

    printf("European put: %.4f\n", pricer.price(euro));                 // ~5.57  (Black-Scholes 5.5735)
    printf("American put: %.4f\n", pricer.price(amer));                 // ~6.25  (early-exercise premium)
    printf("Asian put   : %.4f\n", pricer.priceArithmeticAsian(euro));  // ~3.38
    return 0;
}
```

### Compute Greeks

```cpp
#include "greeks/greeks.hpp"

Greeks greeks(euro, pricer);   // any instrument the pricer can price
greeks.setSeed(20260827u);     // optional: makes the run reproducible

GreekSet g = greeks.all();     // or greeks.delta(), greeks.vega(), ...
```

### Link against it

```bash
g++ -std=c++11 -O3 -fopenmp -Iinclude demo.cpp -L. -lIto -o demo
```

---

## Architecture

```
include/
├── core/
│   ├── constants.hpp        OptionConstants: defaults, FD bump sizes, EPSILON
│   └── ticker.hpp           Ticker — the underlying asset: (S0, sigma)
├── instruments/
│   ├── base_option.hpp      BaseOption — the contract: Ticker + K, r, T, isPut, payoff()
│   ├── european_option.hpp  EuropeanOption : BaseOption  (exercise-style tag type)
│   └── american_option.hpp  AmericanOption : BaseOption  (exercise-style tag type)
├── pricers/
│   └── base_pricer.hpp      BasePricer — path generation and all three pricing methods
├── greeks/
│   └── greeks.hpp           Greeks — finite-difference sensitivities
└── math/
    └── lsmc_utils.hpp       MathUtils — Laguerre basis, least squares, Gaussian elimination

src/                         mirrors include/
deprecated/                  pre-refactor pricers, NOT compiled
```

Three concerns are kept apart deliberately:

1. **`Ticker`** — what the *asset* is (spot, volatility). Reusable across contracts.
2. **`BaseOption`** — what the *contract* is (strike, rate, maturity, put/call) and how it pays off. `payoff()` is the only place call/put asymmetry lives, so every pricing method supports both for free.
3. **`BasePricer`** — *how* you value it (paths, time steps). Discretisation belongs to the numerical method, not to the contract.

`EuropeanOption` and `AmericanOption` carry no data of their own. They exist so `BasePricer::price` can overload on exercise style.

---

## Numerics

### Path generation

`generateGBMPaths` uses the exact GBM solution rather than an Euler step:

```
S(t+dt) = S(t) * exp[ (r - 0.5*sigma^2)*dt + sigma*sqrt(dt)*Z ],   Z ~ N(0,1)
```

so there is **no time-discretisation bias for European payoffs** — only sampling error. Path-dependent payoffs (Asian, American) still carry discretisation error from monitoring at `N` discrete dates.

Each OpenMP thread gets its own `std::mt19937`, seeded from a shared draw plus `omp_get_thread_num()`. Generation is therefore thread-safe, and reproducible only when a seed is fixed via `BasePricer::setSeed`.

### Longstaff–Schwartz

Backward induction over the path set:

1. Set terminal cash flows to the terminal payoff.
2. Walk backwards from `N-1` to `1`. At each date, take only the in-the-money paths.
3. Regress discounted future cash flow on `{1, L0(x), L1(x), L2(x)}` where `x = S/K`. Scaling by moneyness keeps the exponential basis numerically well behaved.
4. Exercise where immediate payoff beats the fitted continuation value; otherwise roll the cash flow forward.
5. Discount the surviving cash flows at `t=1` back to today.

The regression is skipped when four or fewer paths are in the money — with four unknowns the normal-equations matrix would be singular.

---

## Greeks

Validated against closed-form Black–Scholes. European put, S=K=100, r=5%, σ=20%, T=1, 200k paths:

| | `Greeks` | Black–Scholes | error |
| --- | ---: | ---: | ---: |
| price | 5.60079 | 5.57353 | +0.49% |
| delta | −0.36480 | −0.36317 | −0.45% |
| gamma | 0.01874 | 0.01876 | −0.11% |
| vega | 37.68755 | 37.52403 | +0.44% |
| rho | −42.07588 | −41.89046 | −0.44% |
| theta | −1.66527 | −1.65788 | −0.45% |

### Common random numbers are not optional

Each Greek is a difference of two Monte Carlo estimates divided by a small bump. If the two legs draw independent randoms, their sampling error does not cancel — and dividing by the bump amplifies whatever is left. Measured here, delta with an absolute 0.01 bump:

| | mean | std dev | observed range |
| --- | ---: | ---: | --- |
| Independent draws | −0.207 | 2.027 | [−4.09, +2.29] |
| Common random numbers | −0.3631 | 0.0014 | [−0.365, −0.361] |

A put delta is confined to [−1, 0]; the first row gets the *sign* wrong about a third of the time. `Greeks` therefore pins the seed across both legs by default. `setCommonRandomNumbers(false)` disables it, which is useful only for reproducing the table above. A `SeedGuard` restores the pricer's seed state on scope exit, so computing Greeks never leaves a caller's pricer reconfigured.

### Bump sizes and units

Spot bumps are **relative** (a fraction of S0); the rest are absolute. Gamma takes a larger default bump than delta because the second difference divides by `h²`, amplifying noise quadratically — at 1% of spot, gamma carries ~16% noise even with CRN, versus ~0.7% at 5%.

| Greek | default bump | units of the result |
| --- | --- | --- |
| delta | 1% of spot | per unit of spot |
| gamma | 5% of spot | per unit of spot squared |
| vega | 0.01 vol | per 1.0 of vol (÷100 for per vol point) |
| rho | 0.001 | per 1.0 of rate (÷10000 for per bp) |
| theta | 1/365 year | per year (÷365 for per day) |

Override with `GreekBumps` via the constructor's third argument or `setBumps`. Vega and theta fall back to a one-sided difference when a bump would push σ or T through zero.

`all()` costs 10 pricer calls; individual accessors cost 1–3 each.

### Why only the constructor is a template

`BasePricer::price` is a set of compile-time **overloads** — the compiler picks one from the *written* type of the argument. A `BaseOption` matches neither candidate, so the obvious design does not compile:

```
error: no matching function for call to 'BasePricer::price(BaseOption&)'
note: candidate 1: 'double BasePricer::price(const EuropeanOption&)'
note:   no known conversion for argument 1 from 'BaseOption' to 'const EuropeanOption&'
```

The type is known in exactly one place: the construction site. So the constructor — and only the constructor — is templated. It captures the concrete instrument in a `std::function<double(const BaseOption&)>` that transplants bumped terms onto a correctly-typed copy and calls `price()`. That single function absorbs the entire overload problem.

Everything else is ordinary code in `src/greeks/greeks.cpp`. The contract lives in a plain `BaseOption option_` member, and bumped scenarios are built with `BaseOption` and `Ticker` themselves — `setTicker`, `setK`, `setR`, `setT`, and their validation. Nothing about the contract is re-modelled.

Two payoffs: `Greeks` is one plain type, so `std::vector<Greeks>` holds European and American positions side by side; and because the lambda copies the *concrete* option each call, an instrument that later carries its own state (a barrier level, an averaging schedule) keeps it, where a sliced `BaseOption` would not. The cost is one indirect call per pricing — nothing against a 200k-path simulation.

---

## Build system

Two independent builds producing the same artifact. **CMake is the primary one**; the Makefile is a fallback.

### CMakeLists.txt

| Directive | Why it is there |
| --- | --- |
| `cmake_minimum_required(VERSION 3.10)` | Floor for the `OpenMP::OpenMP_CXX` imported target. |
| `set(CMAKE_CXX_STANDARD 11)` + `..._REQUIRED ON` | C++11, and *fail* rather than silently downgrade. |
| `find_package(OpenMP REQUIRED)` | Locates the compiler's OpenMP flags and runtime. |
| `if(MSVC) ... else()` | MSVC and GCC/Clang spell flags completely differently (`/O2` vs `-O3`). |
| `add_compile_options(/wd4267)` | Silences MSVC's `size_t` → `int` narrowing warning. |
| `file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")` | Auto-discovers sources. `CONFIGURE_DEPENDS` re-checks the glob every build, so new files appear without re-running configure. It globs only `src/`, which is why `deprecated/` is excluded for free. |
| `add_library(Ito STATIC ${SOURCES})` | A static archive: not an executable, not a DLL. |
| `target_link_libraries(Ito PUBLIC OpenMP::OpenMP_CXX)` | `PUBLIC` propagates OpenMP to consumers — they never type `-fopenmp`. |
| `target_include_directories(Ito PUBLIC include)` | Same for headers: consumers get `include/` automatically. |

**Multi-config caveat.** The Visual Studio generator is *multi-config*: `CMAKE_BUILD_TYPE` is ignored at configure time and the configuration is chosen at build time with `--config Release`. Makefiles and Ninja are the reverse — `-DCMAKE_BUILD_TYPE=Release` at configure time, no `--config`. Getting this backwards is the usual cause of an accidentally unoptimised build.

### Makefile

```make
SOURCES = $(shell find $(SRCDIR) -name '*.cpp')
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
```

`find` collects sources; `patsubst` rewrites `src/pricers/base_pricer.cpp` into `obj/pricers/base_pricer.o`, mirroring the tree under `obj/`.

```make
$(TARGET): $(OBJECTS)
	ar rcs $@ $^
```

`ar rcs` builds the archive — `r` insert, `c` create quietly, `s` write the symbol index. `$@` is the target, `$^` all prerequisites.

```make
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
```

A pattern rule: any `.o` is built from the matching `.cpp`. `mkdir -p $(dir $@)` creates the nested output directory, since `patsubst` invents subdirectories that do not exist yet. `$<` is the first prerequisite, `$@` the output. Make recompiles only when the `.cpp` is newer — that is the incremental build.

Two limitations worth knowing:

- **`make debug` straight after `make` does nothing.** The target-specific `CXXFLAGS += -g` applies only to the `debug` target, but the objects are already newer than their sources, so Make skips them. Run `make clean debug`.
- **No header dependency tracking.** Make only knows a `.o` depends on its `.cpp`. Edit a header and nothing rebuilds — run `make clean`. CMake tracks this correctly.

The Makefile also needs a POSIX shell (`find`, `mkdir -p`), so Git Bash, MSYS2, WSL, Linux or macOS. Under `cmd.exe`, Windows' own `find.exe` is a text search tool and the build silently produces nothing. Use CMake on native Windows.

Note the Makefile compiles with `-Wall -Wextra` while the CMake build does not — the Makefile is the stricter of the two.

---

## Known gaps

Roughly in priority order.

1. **`AmericanOption` carries an unused RNG.** `gen` and `norm_dist` are leftovers from when the instrument priced itself; path generation now lives in `BasePricer::generateGBMPaths`. Marked `TODO(unused)` in the header. Removing them also lets `<random>` drop out of that header.
2. **LSMC has a quadratic hot spot.** The "roll forward the out-of-the-money paths" step linearly scans `inMoneyIndices` for every path (`src/pricers/base_pricer.cpp`), giving O(paths × in-the-money paths) per time step. A `std::vector<char>` membership mask makes it O(paths). Measured: at 100k paths an American price drops from 22.0s to 1.5s, and scaling goes from ~4.7× per doubling to ~2.2×. This currently caps usable path counts, and matters more now that a full `all()` costs ten pricings.
3. **`BasePricer::price` overloads are not virtual**, despite the `virtual ~BasePricer()`. The polymorphism the destructor implies is not wired up. Decide: plain overloads on a concrete class, or a virtual `price(const BaseOption&)` interface — the latter would also let `Greeks` drop its templated constructor entirely.
4. **`/openmp` is set twice on MSVC** — once in `CMAKE_CXX_FLAGS_RELEASE`, once via `OpenMP::OpenMP_CXX`. Harmless, but drop it from the manual flags.
5. **No tests, no benchmarks, no CI.** The Black–Scholes comparison above is the obvious first test to commit.
6. **Gamma bump sizing** may need tuning away from the 5%-of-spot default for short-dated or far out-of-the-money contracts; the bias/variance tradeoff moves with the contract.
7. **Everything is in the global namespace.** `Ticker` in particular will collide. `namespace ito { ... }` is the natural fix.
8. **Only first- and second-order Greeks.** No vanna, volga, charm, or speed.

---

## Renaming

The name appears in exactly five places, all in the two build files — nothing in `include/` or `src/` hard-codes it:

| File | Line |
| --- | --- |
| `CMakeLists.txt` | 2, 29, 33, 36 |
| `Makefile` | 5 |

```bash
sed -i 's/Ito/NewName/g' CMakeLists.txt Makefile && rm -rf build
```

Delete `build/` afterwards — a stale `CMakeCache.txt` pins the old project name and source directory.

---

## Requirements

- A C++11 compiler with OpenMP: MSVC 2019+, GCC 7+, or Clang with `libomp`
- CMake 3.10+ *(or)* GNU Make with a POSIX shell
- No third-party libraries
