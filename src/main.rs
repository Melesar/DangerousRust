#![allow(
    non_upper_case_globals,
    non_camel_case_types,
    non_snake_case
)]

use std::arch::x86_64::*;
use std::f64::consts::PI;


const STARTING_CONDITION : [body; BODIES_COUNT] = [
    body {    // Sun
        mass: SOLAR_MASS,
        position: [0.; 3],
        velocity: [0.; 3],
    },
    body {    // Jupiter
        position: [
            4.84143144246472090e+00,
            -1.16032004402742839e+00,
            -1.03622044471123109e-01
        ],
        velocity: [
            1.66007664274403694e-03 * DAYS_PER_YEAR,
            7.69901118419740425e-03 * DAYS_PER_YEAR,
            -6.90460016972063023e-05 * DAYS_PER_YEAR
        ],
        mass: 9.54791938424326609e-04 * SOLAR_MASS,
    },
    body {    // Saturn
        position: [
            8.34336671824457987e+00,
            4.12479856412430479e+00,
            -4.03523417114321381e-01
        ],
        velocity: [
            -2.76742510726862411e-03 * DAYS_PER_YEAR,
            4.99852801234917238e-03 * DAYS_PER_YEAR,
            2.30417297573763929e-05 * DAYS_PER_YEAR
        ],
        mass: 2.85885980666130812e-04 * SOLAR_MASS
    },
    body {    // Uranus
        position: [
            1.28943695621391310e+01,
            -1.51111514016986312e+01,
            -2.23307578892655734e-01
        ],
        velocity: [
            2.96460137564761618e-03 * DAYS_PER_YEAR,
            2.37847173959480950e-03 * DAYS_PER_YEAR,
            -2.96589568540237556e-05 * DAYS_PER_YEAR
        ],
        mass: 4.36624404335156298e-05 * SOLAR_MASS
    },
    body {    // Neptune
        position: [
            1.53796971148509165e+01,
            -2.59193146099879641e+01,
            1.79258772950371181e-01
        ],
        velocity: [
            2.68067772490389322e-03 * DAYS_PER_YEAR,
            1.62824170038242295e-03 * DAYS_PER_YEAR,
            -9.51592254519715870e-05 * DAYS_PER_YEAR
        ],
        mass: 5.15138902046611451e-05 * SOLAR_MASS
    }
];

#[repr(C)]
struct body {
    position: [f64; 3],
    velocity: [f64; 3],
    mass: f64,
}

const SOLAR_MASS : f64 = 4_f64 * PI * PI;
const DAYS_PER_YEAR : f64 = 365.24;
const BODIES_COUNT : usize = 5;
const INTERACTIONS_COUNT : usize = BODIES_COUNT * (BODIES_COUNT - 1) / 2;
const ROUNDED_INTERACTIONS_COUNT : usize = INTERACTIONS_COUNT + INTERACTIONS_COUNT % 2;

#[repr(C)]
#[derive(Copy, Clone)]
union Interactions {
    scalars: [f64; ROUNDED_INTERACTIONS_COUNT],
    vectors: [__m128d; ROUNDED_INTERACTIONS_COUNT / 2],
}

impl Interactions {
    // Safety: the in-memory representation of `f64` and `__m128d` is
    // compatible, so accesses to the union members is safe in any
    // order.
    fn as_scalars(&mut self) -> &mut [f64; ROUNDED_INTERACTIONS_COUNT] {
        unsafe {
            &mut self.scalars
        }
    }

    // Safety: the in-memory representation of `f64` and `__m128d` is
    // compatible, so accesses to the union members is safe in any
    // order.
    fn as_vectors(&mut self) -> &mut [__m128d; ROUNDED_INTERACTIONS_COUNT / 2] {
        unsafe {
            &mut self.vectors
        }
    }
}

fn offset_Momentum(bodies: &mut [body; BODIES_COUNT]) {
    for i in  0..BODIES_COUNT {
        for m in 0..3 {
            bodies[0].velocity[m] -= 
                bodies[i].velocity[m] *
                bodies[i].mass / SOLAR_MASS;
        }
    }
}

fn output_Energy(bodies: &[body; BODIES_COUNT]) {
    let mut energy = 0_f64;
    for i in 0..BODIES_COUNT {
        energy += 0.5 * bodies[i].mass * (
            bodies[i].velocity[0] * bodies[i].velocity[0] +
            bodies[i].velocity[1] * bodies[i].velocity[1] +
            bodies[i].velocity[2] * bodies[i].velocity[2]
        );

        for j in i + 1..BODIES_COUNT {
            let mut position_Delta = [0.; 3];
            for m in 0..3 {
                position_Delta[m] = bodies[i].position[m] - bodies[j].position[m];
            }

            energy -= bodies[i].mass * bodies[j].mass / f64::sqrt(
                position_Delta[0] * position_Delta[0] +
                position_Delta[1] * position_Delta[1] +
                position_Delta[2] * position_Delta[2]
            );
        }
    }

    println!("{:.9}", energy);
}

#[cfg(target_feature = "sse2")]
fn advance(bodies: &mut [body; BODIES_COUNT],
                  position_Deltas: &mut [Interactions; 3],
                  magnitudes: &mut Interactions) {

    {
        let mut k : usize = 0;
        for i in 0..BODIES_COUNT - 1 {
            for j in i + 1..BODIES_COUNT {
                for m in 0..3 {
                    position_Deltas[m].as_scalars()[k] = bodies[i].position[m] - 
                        bodies[j].position[m];
                }

                k += 1;
            }
        }
    }

    for i in 0..ROUNDED_INTERACTIONS_COUNT/2 {
        // Safety: this code is only compiled for processors that
        // support SSE2, so the SIMD operations used here are
        // safe.
        let mut position_Delta = unsafe { [_mm_set1_pd(0.); 3] };
        for m in 0..3 {
            position_Delta[m] = position_Deltas[m].as_vectors()[i];
        }

        // Safety: this code is only compiled for processors that
        // support SSE2, so the SIMD operations used here are
        // safe.
        let distance_Squared : __m128d = unsafe {
            _mm_add_pd(
                _mm_mul_pd(position_Delta[0], position_Delta[0]),
                _mm_add_pd(
                    _mm_mul_pd(position_Delta[1], position_Delta[1]),
                    _mm_mul_pd(position_Delta[2], position_Delta[2]),
                    )
            )
        };

        // Safety: this code is only compiled for processors that
        // support SSE2, so the SIMD operations used here are
        // safe.
        let mut distance_Reciprocal : __m128d = unsafe {
            _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(distance_Squared)))
        };

        for _ in 0..2 {
            // Safety: this code is only compiled for processors that
            // support SSE2, so the SIMD operations used here are
            // safe.
            distance_Reciprocal = unsafe {
                _mm_sub_pd(
                    _mm_mul_pd(distance_Reciprocal, _mm_set1_pd(1.5)),
                    _mm_mul_pd(
                        _mm_mul_pd(
                            _mm_mul_pd(
                                _mm_set1_pd(0.5),
                                distance_Squared
                            ),
                            distance_Reciprocal
                        ),
                        _mm_mul_pd(distance_Reciprocal, distance_Reciprocal)
                      )
                )
            };
        }

        // Safety: this code is only compiled for processors that
        // support SSE2, so the SIMD operations used here are
        // safe.
        magnitudes.as_vectors()[i] = unsafe {
            _mm_mul_pd(
                _mm_div_pd(
                    _mm_set1_pd(0.01),
                    distance_Squared
                ),
                distance_Reciprocal
            )
        };
    }

    {
        let mut k : usize = 0;
        for i in 0..BODIES_COUNT {
            for j in i + 1..BODIES_COUNT {
                let i_mass_magnitude = bodies[i].mass * magnitudes.as_scalars()[k];
                let j_mass_magnitude = bodies[j].mass * magnitudes.as_scalars()[k];
                for m in 0..3 {
                    bodies[i].velocity[m] -= position_Deltas[m].as_scalars()[k] * j_mass_magnitude;
                    bodies[j].velocity[m] += position_Deltas[m].as_scalars()[k] * i_mass_magnitude;
                }
                k += 1;
            }
        }
    }

    for i in 0..BODIES_COUNT {
        for m in 0..3 {
            bodies[i].position[m] += 0.01 * bodies[i].velocity[m];
        }
    }
}


fn main() {
    let mut solar_Bodies = STARTING_CONDITION;
    let mut position_Deltas = [Interactions {scalars: [0.; ROUNDED_INTERACTIONS_COUNT]}; 3];
    let mut magnitudes = Interactions {scalars: [0.; ROUNDED_INTERACTIONS_COUNT]};

    offset_Momentum(&mut solar_Bodies);
    output_Energy(&solar_Bodies);

    let c = std::env::args().nth(1).unwrap().parse().unwrap();
    for _ in 0..c {
        advance(&mut solar_Bodies, &mut position_Deltas, &mut magnitudes);
    }
    
    output_Energy(&solar_Bodies);
}
