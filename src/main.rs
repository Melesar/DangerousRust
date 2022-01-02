use std::f64::consts::PI;

const STARTING_CONDITION : [Body; BODIES_COUNT] = [
    Body {    // Sun
        mass: SOLAR_MASS,
        position: [0.; 3],
        velocity: [0.; 3],
    },
    Body {    // Jupiter
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
    Body {    // Saturn
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
    Body {    // Uranus
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
    Body {    // Neptune
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

#[derive(Clone, Debug)]
struct Body {
    position: [f64; 3],
    velocity: [f64; 3],
    mass: f64,
}

const SOLAR_MASS : f64 = 4_f64 * PI * PI;
const DAYS_PER_YEAR : f64 = 365.24;
const BODIES_COUNT : usize = 5;
const INTERACTIONS_COUNT : usize = BODIES_COUNT * (BODIES_COUNT - 1) / 2;

fn sqr(x: f64) -> f64 {
    x * x
}

fn offset_momentum(bodies: &mut [Body; BODIES_COUNT]) {
    let (sun, planets) = bodies.split_first_mut().unwrap();
    sun.velocity = [0.; 3];

    for planet in planets {
        for m in 0..3 {
            sun.velocity[m] -= planet.velocity[m] * planet.mass / SOLAR_MASS;
        }
    }
}

fn output_energy(bodies: &[Body; BODIES_COUNT]) {
    let mut energy = 0_f64;
    for (i, body) in bodies.iter().enumerate() {
        energy += 0.5 * body.mass * (
            sqr(body.velocity[0]) +
            sqr(body.velocity[1]) +
            sqr(body.velocity[2])
        );

        for body2 in &bodies[i + 1..BODIES_COUNT] {
            energy -= body.mass * body2.mass / f64::sqrt(
                sqr(body.position[0] - body2.position[0]) +
                sqr(body.position[1] - body2.position[1]) +
                sqr(body.position[2] - body2.position[2])
            );
        }
    }

    println!("{:.9}", energy);
}

fn advance(bodies: &mut [Body; BODIES_COUNT]) {

    let mut position_deltas = [[0.; 3]; INTERACTIONS_COUNT];
    {
        let mut k : usize = 0;
        for i in 0..BODIES_COUNT - 1 {
            for j in i + 1..BODIES_COUNT {
                for m in 0..3 {
                    position_deltas[k][m] = bodies[i].position[m] - 
                        bodies[j].position[m];
                }

                k += 1;
            }
        }
    }

    let magnitudes = {
        let mut magnitudes = [0.; INTERACTIONS_COUNT];
        for (i, pd) in position_deltas.iter().enumerate() {
            let distance_squared = sqr(pd[0]) + sqr(pd[1]) + sqr(pd[2]);
            magnitudes[i] = 0.01 / (distance_squared * distance_squared.sqrt());
        }
        magnitudes
    };

    {
        let mut k : usize = 0;
        for i in 0..BODIES_COUNT {
            for j in i + 1..BODIES_COUNT {
                let i_mass_magnitude = bodies[i].mass * magnitudes[k];
                let j_mass_magnitude = bodies[j].mass * magnitudes[k];
                for (m, pd) in position_deltas[k].iter().enumerate() {
                    bodies[i].velocity[m] -= pd * j_mass_magnitude;
                    bodies[j].velocity[m] += pd * i_mass_magnitude;
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
    let mut solar_bodies = STARTING_CONDITION;

    offset_momentum(&mut solar_bodies);
    output_energy(&solar_bodies);

    let c = std::env::args().nth(1).unwrap().parse().unwrap();
    for _ in 0..c {
        advance(&mut solar_bodies);
    }
    
    output_energy(&solar_bodies);
}
