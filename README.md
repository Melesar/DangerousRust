# Dangerous Rust

This is my walkthrough the excellent project of [Cliff L. Biffle](https://cliffle.com/about/) - [Learn Rust the Dangerous Way](https://cliffle.com/p/dangerust/). 

Here I try to follow along the tutorial, learning new techniques and running the benchmarks myself. Each step of the tutorial is represented by a commit, so you can how the program evolved looking at the commit log. The final, [bonus](https://cliffle.com/p/dangerust/6/) step is implemented on a separate branch: _rusty_

_**Spoiler alert**_: Running the benchmarks on the 5-th and the final version of the program against `gcc` gave the following results:

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `./nbody.gcc-8 50000000` | 2.380 ± 0.028 | 2.352 | 2.425 | 1.00 |
| `./nbody-5.bench 50000000` | 2.573 ± 0.005 | 2.567 | 2.583 | 1.08 ± 0.01 |
| `./nbody-6.bench 50000000` | 2.586 ± 0.028 | 2.569 | 2.645 | 1.09 ± 0.02 |

To me this is surprising, since it's quite different from the author's results. First thing is that `gcc` managed to slightly outperform Rust in both cases. And the second is that the final version didn't give any significant advantage over the 5-th, which was the case for the author.

My version of gcc: `gcc (Ubuntu 11.2.0-7ubuntu2) 11.2.0`
My version of rustc: `rustc 1.56.1 (59eed8a2a 2021-11-01)`

