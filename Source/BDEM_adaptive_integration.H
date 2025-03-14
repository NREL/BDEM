
// Step 1: Compute forces at current time
BL_PROFILE_VAR( "FORCE_CALC", forceCalc );
{

    bpc.computeForces(
        dt,
        EBtools::ebfactory,
        EBtools::lsphi,
        specs.do_heat_transfer,
        specs.walltemp_vardir,
        specs.walltemp_polynomial.data(),
        EBtools::ls_refinement,
        specs.stl_geom_present,
        specs.contact_law,
        steps,
        specs.gravity,
        specs.stls,
        time,
        specs.bonded_sphere_particles,
        specs.liquid_bridging,
        specs.particle_cohesion,
        specs.init_force,
        specs.init_force_dir,
        specs.init_force_comp,
        cb_force,
        cb_torq,
        specs.cb_dir,
        specs.drag_model,
        specs.solve_fibrillation,
        specs.fib_alpha,
        specs.fib_beta,
        specs.fib_ced
    );
}
BL_PROFILE_VAR_STOP( forceCalc );

// Step 2: Advance velocity and position with a low order scheme
BL_PROFILE_VAR( "MOVE_PART", movepart );
specs.verlet_scheme = 1;
bpc.moveParticles( dt, specs.do_chemistry, specs.minradius_frac, specs.verlet_scheme, 1 );
BL_PROFILE_VAR_STOP( movepart );


// Step 3: Recompute forces
BL_PROFILE_VAR( "FORCE_CALC", forceCalc );
{

    bpc.computeForces(
        dt,
        EBtools::ebfactory,
        EBtools::lsphi,
        specs.do_heat_transfer,
        specs.walltemp_vardir,
        specs.walltemp_polynomial.data(),
        EBtools::ls_refinement,
        specs.stl_geom_present,
        specs.contact_law,
        steps,
        specs.gravity,
        specs.stls,
        time,
        specs.bonded_sphere_particles,
        specs.liquid_bridging,
        specs.particle_cohesion,
        specs.init_force,
        specs.init_force_dir,
        specs.init_force_comp,
        cb_force,
        cb_torq,
        specs.cb_dir,
        specs.drag_model,
        specs.solve_fibrillation,
        specs.fib_alpha,
        specs.fib_beta,
        specs.fib_ced
    );
}
BL_PROFILE_VAR_STOP( forceCalc );

// Step 4: Apply high order correction to velocity and position
BL_PROFILE_VAR( "MOVE_PART", movepart );
bpc.moveParticles( dt, specs.do_chemistry, specs.minradius_frac, specs.verlet_scheme, 2 );
BL_PROFILE_VAR_STOP( movepart );

