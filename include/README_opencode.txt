What processes the .fipps files
- source/input_tf.f90 (2223 lines) — reads control.fipps, then, depending on the keywords set there, reads ~45 sibling .fipps files (nodes, loads, elements, materials, failure criteria, properties, boundary conditions, laminates, couplings, contact, aero-coupling) into the fe_simulation type fesim.
- source/input_process_line.f90 — the keyword/value parser for control.fipps.

The API: include/fipps_api.f90 (new module fipps_api)
Replaces every *.fipps file with a typed subroutine taking arrays as input arguments — no file I/O, no goto, no read statements. Design:
- One subroutine per file — each takes fesim (intent(inout)), one dummy array per file column, and an err output. E.g. fipps_api_nodes(fesim, cids, coords, err) replaces nodes.fipps.
- Each call sets the matching is_<section> flag, so the control-file keyword mechanism is fully replaced by fipps_api_control(...) (all options are optional named arguments) plus the section calls themselves.
- fipps_api_check(fesim, err) reproduces the control.fipps consistency checks (input_tf.f90:143-174).
- Faithful to the original semantics: nid/eid numbering, eid running across element types (beam2→quad8→lsolid20), sfac=1 default for loads, mpc_type=0, lay=0 for pcomp, dfac clamping, fid<1→0, angle deg→rad conversion, coord column normalization, contact rows grouped per element.
- Error codes: 0 ok, 1 bad args, 2 bad angle unit, 3 section already provided, 4 precondition missing.
- failure.fipps (mixed type-names per row) → 11 dedicated subroutines fipps_api_failure_{tresca,mises,maxpstress,puck,hill,norris,fibre,maxstrain,cuntze,maxstrain3d,tsaiwu3d}.
- Jagged structures handled: MPC master DOFs via a (nmax,n) matrix + per-row counts; lam20 per-subcase material IDs via mat20ids(nsub,n); the two aero-coupling files per dimension each collapsed into one call.
- Added fipps_api.o to the module list in include/Makefile.

Verification
- Compiles cleanly with ifx (Intel oneAPI 2025.1.1) alongside all existing modules.
- A runtime smoke-test driver built a mixed beam2+quad8+lsolid20 model purely through the API and verified data landed correctly in fesim (node/element IDs, angle conversion, MPC grouping, contact grouping, flags, err codes) — all checks passed.

Usage flow replacing input_tf:
call init_default(fesim)
call fipps_api_control(fesim, sol=1, outputVTK=.true., err=err)
call fipps_api_nodes   (fesim, cids, coords, err)
! ... one call per section ...
call fipps_api_check   (fesim, err)
call init_values(fesim)

One constraint to keep in mind: element sections must be called in the order beam2 → quad8 → lsolid20 (same as the original file reader) so the global element IDs stay identical.
