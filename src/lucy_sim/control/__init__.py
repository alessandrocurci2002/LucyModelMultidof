"""Controllers in the master: LQG-PI of Gawronski (:mod:`.lqg`). The cascade of the drive is in :mod:`lucy_sim.plant`.

A controller for the simulators needs ``initial_state()``, ``output(state, r)`` (velocity command for the reference
``r``) and ``update(state, y, u_applied)``, plus ``measure(state, y)`` if it reads the measured position before
computing the command (:class:`lqg.LQGPIMeasured`)."""
