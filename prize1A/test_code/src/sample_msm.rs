




use std::io::Write;
use std::{fs::File, mem::transmute};

use ark_ec::msm::VariableBaseMSM;
use ark_ec::short_weierstrass_jacobian::{GroupAffine, GroupProjective};
use ark_ec::{ProjectiveCurve, SWModelParameters};
use ark_ff::UniformRand;
use ark_ff::{Field, PrimeField, SquareRootField};
use rand_chacha::{rand_core::SeedableRng, ChaCha20Rng};

use crate::Args;

#[derive(Debug)]
struct EdwardsPoint<F> {
    x: F,
    y: F,
    t: F,
}

impl<F> EdwardsPoint<F> {
    fn new(x: F, y: F, t: F) -> Self {
        Self { x, y, t }
    }
}

fn weierstrass_to_edwards<F: PrimeField + SquareRootField>(x: &F, y: &F) -> EdwardsPoint<F> {
    // Constants for the conversion
    let a = F::zero(); // Weierstrass a parameter
    let b = F::one();  // Weierstrass b parameter
    
    // Montgomery conversion parameters
    let alpha = -F::one();
    let s = (F::from(3u64) * alpha.square() + a).sqrt().unwrap().inverse().unwrap();
    let a = F::from(3u64) * alpha * s;
    let b = s;
    
    // Convert to Montgomery form
    let mx = (*x) * b - a / F::from(3u64);
    let my = (*y) * b;
    
    // Convert Montgomery to Edwards
    let ex = mx / my;
    let ey = (mx - F::one()) / (mx + F::one());
    let et = ex * ey; // t-coordinate for extended Edwards form
    
    EdwardsPoint::new(ex, ey, et)
}

pub(crate) fn sample_msm<P: SWModelParameters<BaseField = F>, F: PrimeField + SquareRootField>(args: &Args) {
    let mut rng = ChaCha20Rng::seed_from_u64(args.seed);

    let scalars: Vec<P::ScalarField> = (0..args.degree)
        .map(|_| P::ScalarField::rand(&mut rng))
        .collect::<Vec<_>>();
    let scalars_bigint = scalars.iter().map(|x| x.into_repr()).collect::<Vec<_>>();

    let bases: Vec<GroupAffine<P>> = (0..args.degree)
        .map(|_| GroupProjective::<P>::rand(&mut rng).into_affine())
        .collect::<Vec<_>>();

    let res = VariableBaseMSM::multi_scalar_mul(&bases, &scalars_bigint).into_affine();

    // for (i, e) in scalars.iter().enumerate() {
    //     println!("{} {:?}", i, e)
    // }

    // for (i, e) in bases.iter().enumerate() {
    //     println!("{} {:?}", i, e)
    // }

    // scalars
    {
        let filename = format!(
            "{}_scalars_dim_{}_seed_{}.csv",
            F::size_in_bits(),
            scalars.len(),
            args.seed,
        );
        let mut file = File::create(filename).unwrap();

        for (i, s) in scalars.iter().enumerate() {
            file.write_all(format!("{}, {}\n", i, printer_256(s)).as_ref())
                .unwrap();
        }
    }

    // base
    {
        let filename = format!(
            "{}_bases_dim_{}_seed_{}.csv",
            F::size_in_bits(),
            scalars.len(),
            args.seed,
        );
        let mut file = File::create(filename).unwrap();

        for (i, s) in bases.iter().enumerate() {
            // Convert each base point to Edwards form
            let edwards_point = weierstrass_to_edwards(&s.x, &s.y);
            
            file.write_all(
                format!(
                    "{}, x, {}, y, {}, t, {}\n",
                    i,
                    printer_384(&edwards_point.x),
                    printer_384(&edwards_point.y),
                    printer_384(&edwards_point.t)
                )
                .as_ref(),
            ).unwrap();
        }
    }

    // result
    {
        let filename = format!(
            "{}_res_dim_{}_seed_{}.csv",
            F::size_in_bits(),
            scalars.len(),
            args.seed,
        );
        let mut file = File::create(filename).unwrap();

        // Convert result to Edwards form
        let edwards_res = weierstrass_to_edwards(&res.x, &res.y);
        
        file.write_all(
            format!(
                "x, {}, y, {}, t, {}\n", 
                printer_384(&edwards_res.x), 
                printer_384(&edwards_res.y),
                printer_384(&edwards_res.t)
            ).as_ref(),
        ).unwrap();
    }
}

fn printer_256<F: Field>(f: &F) -> String {
    let mont = unsafe { transmute::<_, &[u64; 4]>(f) };
    format!(
        "{:0>16x?}, {:0>16x?}, {:0>16x?}, {:0>16x?}, zz: {}",
        mont[3], mont[2], mont[1], mont[0], f
    )
}

fn printer_384<F: Field>(f: &F) -> String {
    let mont = unsafe { transmute::<_, &[u64; 6]>(f) };
    format!(
        "{:0>16x?}, {:0>16x?}, {:0>16x?}, {:0>16x?}, {:0>16x?}, {:0>16x?}, zz: {}",
        mont[5], mont[4], mont[3], mont[2], mont[1], mont[0], f
    )
}
