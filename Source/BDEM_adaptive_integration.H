// Step 0: Compute forces if first step
//! Only for first step!
if ( steps == 0 )
{
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
}


bool redo_step = false;

do
{
    // Step 1: Advance velocity and position with a low order scheme
    BL_PROFILE_VAR( "MOVE_PART", movepart );
    bpc.moveParticles( dt, specs.do_chemistry, specs.minradius_frac, 1, 1 );
    BL_PROFILE_VAR_STOP( movepart );


    // Step 2: Recompute forces
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
            0., // Do not evolve the equation for the fraction of fibrils (only update the force)
            0., // Do not evolve the equation for the fraction of fibrils (only update the force)
            specs.fib_ced
        );
    }
    BL_PROFILE_VAR_STOP( forceCalc );

    // Step 3: Apply high order correction to velocity and position
    BL_PROFILE_VAR( "MOVE_PART", movepart );
    bpc.moveParticles( dt, specs.do_chemistry, specs.minradius_frac, 1, 2 );
    BL_PROFILE_VAR_STOP( movepart );

    redo_step = bpc.updateTimeStep(dt, steps, specs.dt_min, specs.dt_max, specs.dt_error_jump_ratio, specs.dt_error_target );

} while (redo_step);



// Step 4: Compute updated forces
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

