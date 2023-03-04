The Pollard Rho algorithm is the simplest factorization algorithm that's better than trial division in many practical scenarios.
It is a randomized algorithm with a heuristic average runtime of `O(sqrt(p))` where `p` is the smallest prime factor of `n`.

However, this average runtime isn't exactly correct, and the full story is pretty interesting.

Let's look at the assumptions that lead to this runtime, the cases where these assumptions are wrong, and some more exact numbers.

There are three key facts/assumptions that underpin Pollard Rho:

1. Any sequence of `n + 1` numbers mod `n` must contain a repeat, aka the pidgeonhole principle.  There are a couple of refinements to this statement:
   - for recursively defined sequences that depend only on the previous term, this means the entire sequence must repeat with a period of at most `n`
   - for random sequences, the expected number of random numbers mod `n` that can be generated before seeing a repeat is about `sqrt(ln(2)n)`.
     This comes from analysis of the so-called birthday paradox, see https://en.wikipedia.org/wiki/Birthday_problem#Approximations
2. For a good choice of function `g(x)` the corresponding sequence `x_i = g(x_(i-1)) mod n` is pseudorandom.
3. The expected number of steps before the sequence will repeat mod `p|n` is only `sqrt(ln(2)p)`, so with a good cycle detection algorithm we expect to
     find a repeat mod `p` faster than mod `n`

How is it possible to detect a repeat mod `p` when we don't know `p`?
Using a common trick in factoring algorithms: consider `gcd(x_i - x_j mod n, n)`.
If `x_i = x_j mod p`, then `p` will divide this gcd, but otherwise, it will not and the gcd could be 1 or a multiple of other prime factors of `n` if any.

Pollard Rho is usually implemented using Floyd's cycle detection algorithm, which keeps only 2 numbers, `x` and `y`.
Both are started at the same random value, then `x` is advanced one step at a time, ie `x = g(x) mod n`, and `y` is advanced two steps at a time, ie `y = g(g(y)) mod n`.
When `gcd(abs(x - y), n) > 1`, we have found a repeat mod `n`, and, if the gcd is not `n`, mod a factor of `n`, in which case we also have a nontrivial factor of `n`.

However, it's important to notice that Floyd cycle detection doesn't always find a cycle on the first repeated element.
This is a tradeoff: to detect a cycle immediately, we would need to store every single element, but Floyd cycle detection only needs to store 2 elements.
This tradeoff is absolutely worth it, but it does cause Pollard Rho to fail to factor 4 and 25.
In practice this doesn't matter because removing small factors like 2 and 5 with trial division is common practice before moving on to asymptotically better algorithms.
More on that when we look at what happens when we apply Pollard Rho to 4.

When choosing `g(x)`, we are ensuring we have a recursively defined sequence that depends only on the previous term, but we also want it to behave as randomly as possible
while still being simple to compute.

This function could be anything, but `g(x) = x^2 + 1` is by far the most common.
The original algorithm used `g(x) = x^2 - 1` and some implementations use `g(x) = x^2 + a` with random `a`.

We will focus on `g(x) = x^2 + 1`, but the analysis is pretty similar for any polynomial.

Notice that `g(x) = x^2 + 1` produces the same value for `x` and `n - x`, so its range of possible outputs only contains at most `(n + 1)/2` distinct values, not `n`.
This is one big way in which this particular `g` fails to be truly random.

First, let's walk through a couple examples.

### Example 1: Factoring 9

If `n = 9`, there are 9 possible starting values of `x` (see 9.svg for a visual overview):
- `0 -> 1 -> (2 -> 5 -> 8)`
- `1 -> (2 -> 5 -> 8)`
- `(2 -> 5 -> 8)`
- `3 -> 1 -> (2 -> 5 -> 8)`
- `4 -> (8 -> 2 -> 5)`
- `(5 -> 8 -> 2)`
- `6 -> 1 -> (2 -> 5 -> 8)`
- `7 -> (5 -> 8 -> 2)`
- `(8 -> 2 -> 5)`

Notice `g(x)` and `g(9-x)` are always the same, as expected.
So regardless of where we start, we'll eventually fall into the `(2 -> 5 -> 8)` cycle.
However, mod 3, the story is different:
- `0 -> 1 -> (2)`
- `1 -> (2)`
- `(2)`

In other words, we always fall into the `(2)` cycle mod 3.

The Floyd cycle detection algorithm doesn't always find cycles immediately though, so let's analyze what happens in each case (does it find a repeat mod 3 before mod 9 or not):

For mod 3:
- 0: 2 steps
- 1: 1 step
- 2: 1 step

For mod 9:
- 0, 3, and 6: 3 steps
- 1, 4, and 7: 3 steps
- 2, 5, and 8: 3 steps

So regardless of starting point, pollard rho will detect a cycle mod 3 before it detects a cycle mod 9.

The weights for each number of steps before first repeat mod 9 are
- 5 steps with weight 3
- 4 steps with weight 3
- 3 steps with weight 3.

We could compare these with the expected weights for a truly random birthday paradox scenario:
- 2 steps with weight `1/9` ~ 0.1111
- 3 steps with weight `8/9 * 2/9` ~ 0.1975
- 4 steps with weight `8/9 * 7/9 * 3/9` ~ 0.2305
- 5 steps with weight `8/9 * 7/9 * 6/9 * 4/9` ~ 0.2048
- 6 steps with weight `8/9 * 7/9 * 6/9 * 5/9 * 5/9` ~ 0.1423
- 7 steps with weight `8/9 * 7/9 * 6/9 * 5/9 * 4/9 * 6/9` ~ 0.0759
- 8 steps with weight `8/9 * 7/9 * 6/9 * 5/9 * 4/9 * 3/9 * 7/9` ~ 0.0295
- 9 steps with weight `8/9 * 7/9 * 6/9 * 5/9 * 4/9 * 3/9 * 2/9 * 8/9` ~ 0.0075
- 10 steps with weight `8/9 * 7/9 * 6/9 * 5/9 * 4/9 * 3/9 * 2/9 * 1/9` ~ 0.0001

Hopefully this makes the difference clear between
- expected steps until first repeat in a truly random scenario
- expected steps until first repeat in a pollard rho scenario where we us a pseudorandom sequence
- expected steps until floyd cycle detection actually notices the cycle.

### Example 2: Failing to factor 4

Mod 4, there are 4 possible starting values:
- `0 -> (1 -> 2)`
- `(1 -> 2)`
- `(2 -> 1)`
- `3 -> (2 -> 1)`

In the case of 0 and 3, there are 3 steps until the first repeat and 2 steps until the first detection.
In the case of 1 and 2, there are 2 steps until both the first repeat and detection both.

Mod 2, there are 2 possible starting values:
- `(0 -> 1)`
- `(1 -> 0)`

Which, by similarity to 1 and 2 mod 4, we can immediatly tell both repeat and are detected after 2 steps.

So regardless of the starting value, the number of steps to detect a cycle mod 2 is never less than the number of steps to detect a cycle mod 4.

Thus Pollard Rho cannot factor 4, at least with the standard `g(x)` and Floyd (or Brent) cycle detection.

### Example 3: Failing to factor 25

For brevity, we will only look at 0 to 12, since we know that 13-24 are just 25 minus one of those values and thus won't lead to new values of `g(x)`.
We will also skip any number we've seen.
- `0 -> (1 -> 2 -> 5)`
- `3 -> 10 -> (1 -> 2 -> 5)`
- `4 -> 17 -> 15 -> (1 -> 2 -> 5)`
- `6 -> 12 -> 20 -> (1 -> 2 -> 5)`
- `7 -> 0 -> (1 -> 2 -> 5)`
- `8 -> 15 -> (1 -> 2 -> 5)`
- `9 -> 7 -> 0 -> (1 -> 2 -> 5)`
- `11 -> 22 -> 10 -> (1 -> 2 -> 5)`
- `12 -> 20 -> (1 -> 2 -> 5)`

Remember to look at 25.svg for a visual representation!

In this case, we notice that aside from the cycle `(1 -> 2 -> 5)` and the relevant negatives 23 and 24, all of the possible starting points lie on a branch pointing into 1.
There are 4 identical copies of this branch:
- `9/16 -> 7, 7/18 -> 0 -> 1`
- `6/19 -> 12, 12/13 -> 20 -> 1`
- `11/14 -> 22, 22/3 -> 10 -> 1`
- `4/21 -> 17, 17/8 -> 15 -> 1`

Notice how the first number is always `1/4` mod 5, the second is always `2/3` mod 5, the third is `0` mod 5, and the fourth is `1` mod 5.
In other words, each of these branches is isomorphic not only mod 25 but also mod 5, so we just need to analyze one branch.

Also note the graph mod 5:
- `(0 -> 1 -> 2)`

with 3 and 4 also mapping to 0 and 2 respectively, but with nothing pointing to 3 or 4.

Thus,
- if we start at say 9 or 16 mod 25 (or equivalently 6 or 19, 11 or 14, etc), there will be 6 steps until a repeat, and 3 steps until a detection.
   9 and 16 correspond to 4 and 1 mod 5, where there will be 4 or 3 steps until a repeat respectively, and 3 steps until a detection in both cases.
- if we start at 7 or 18 mod 25, there will be 5 steps until a repeat, and 3 steps until a detection.
   7 and 18 correspond to 2 and 3 mod 5, where there will be 3 or 4 steps until a repeat respectively, and 3 steps until a detection in both cases.
- if we start at 0 mod 25, there will be 4 steps until a repeat, and 3 steps until a detection.
   0 corresponds to 0 mod 5, where there will be 3 steps until a repeat, and 3 steps until a detection.
- if we start at 1 mod 25, ther will be 3 steps until both a repeat and a detection, and the same is true for 1 mod 5.
- the reachable subgraphs starting at 2 and 5 are isomorphic to the reachable subgraph starting at 1, so we get 3 steps in all cases mod 25 or 5 again.
- this only leaves 23 and 24, but the reachable subgraphs starting at both are isomorphic to the reachable subgraph starting at 0,
   so there are 4 steps until a repeat mod 25 and 3 steps in the other three cases.

Thus, the number of steps until a repeat is
- 6 and 4 mod 25 and 5 in 4 cases
- 6 and 3 mod 25 and 5 in 4 cases
- 5 and 3 mod 25 and 5 in 4 cases
- 5 and 4 mod 25 and 5 in 4 cases
- 4 and 3 mod 25 and 5 in 4 cases
- 3 and 3 mod 25 and 5 in 3 cases
- 4 and 3 mod 25 and 5 in 2 cases

So if we detected a cycle as soon as it occured, we would find a nontrivial factor in 22/25 cases, but Floyd cycle detection doesn't do that.

The number of steps until a repeat is detected is 3 in all cases, so Pollard Rho never finds a nontrivial factor of 25.

## Analysis of Cycles

2, 3, and 5 give a good idea of how Pollard Rho works, but an inaccurate impression of how the cycles are arranged.
All of them have only one cycle, but most of the time there are many.
However, in the case of 3, the cycle mod 3 is size 1 but mod 9 it is size 3, so keep that interesting fact in mind.

Mod 7 or 49, there are multiple cycles:
- mod 7, we have `(3)` and `(5)`
- mod 49, we have `(3 -> 10)`, `(17 -> 45)`, `(24 -> 38)`, `(31)`, `(19)`, and `(5 -> 26 -> 40 -> 33 -> 12 -> 47)`

The interesting thing is that, aside from 31 and 19, all elements of cycles mod 49 are elements of smaller cycles mod 7.

Notice that, starting at an element of a `k`-cycle, the Floyd algorithm takes `k` steps to detect the cycle.
This is because we can represent the two pointers `x` and `y` in the algorithm as simply numbers mod `k`, where `x = i` and `y = 2i`,
so `x = i = 2i = y mod k <-> i = 0 mod k`.

So whenever we start the Pollard Rho algorithm at a `k`-cyclant (element of a `k`-cycle) mod `p^2` which is a `d`-cyclant mod `p`
where `d < k`, we are guaranteed to find a nontrivial factor.
It's possible to extend this to starting points that reach a cycle in under `d` steps, but it's more important to understand cycles better first.

Although the size of a cycle mod `n` doesn't seem to have many restrictions yet, one key fact is that if `x` is a `k`-cyclant mod `n`, then it will be a `d`-cyclant mod `q`
where `d` divides `k` and `q` divides `n`.

This means that if the largest cycle mod `n` is larger than the largest cycle mod `q`, then starting at any element of the largest cycle mod `n` is guaranteed to find a nontrivial factor of `n`.
Since the maximum cycle size mod `n` is `(n + 1)/2`, this always holds when the largest cycle mod `n` has more than `(q + 1)/2` elements.

Analyzing roots of `g(x) - x`, `g(g(x)) - x`, and further, can give us a little more information.
For example, if `g(x) - x` factors mod `n`, then there are 1-cycles mod `n`, if `g(g(x)) - x` factors mod `n` there are 2-cycles, and so on.

We can also note that, because charmicael's theorem tells us that `x^(lambda + a) = x^a` for any `a > 0`, where `lambda` is the carmichael lambda function,
the sequence of polynomials `x`, `g(x)`, `g(g(x))`, etc, also eventually repeats, although just based on the pidgeonhole principle this could take up to `n^((lambda + 1)/2)`
steps.  Since this represents the number of steps until the tuple of sequence values for each starting point repeats, it should also be bounded above by the maximal
least common multiple of `(n + 1)/2` numbers up to `(n + 1)/2`.  This could be represented as `prod(1 <= j < inf, primorial(((n + 1)/2)^(1/j)))`, and it might even have
a nice representation :)

The maximal cycles mod `p^2` for all primes below 100 are:
- 2: 2
- 3: 3
- 5: 3
- 7: 6
- 11: 20
- 13: 12
- 17: 96
- 19: 9
- 23: 22
- 29: 56
- 31: 30
- 37: 36
- 41: 7

