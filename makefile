# Runs test cases used for paper's results

terminator_toy:
	mpiexec -n 4 python3 test_scripts/terminator_toy.py \
		--ncells_per_edge 16 --dt 450.0 --dirname terminator_toy

galewsky_jet_icosahedron:
	mpiexec -n 120 python3 test_scripts/galewsky_jet.py \
		--ncells_per_edge 64 --dt 150.0 --dumpfreq 3456 --icosahedron

galewsky_jet_cubed_sphere:
	mpiexec -n 120 python3 test_scripts/galewsky_jet.py \
		--ncells_per_edge 100 --dt 150.0 --dumpfreq 3456

baroclinic_wave:
	mpiexec -n 120 python3 test_scripts/baroclinic_wave.py \
		--ncells_per_edge 32

straka_degree_0:
	mpiexec -n 32 python3 test_scripts/straka_degree_0.py \
		--nlayers 128 --dt 0.5 --dumpfreq 450 --dirname straka_degree_0

straka_degree_1:
	mpiexec -n 32 python3 test_scripts/straka_degree_1.py \
		--nlayers 64 --dt 0.5 --dumpfreq 450 --dirname straka_degree_1

moist_3d_bubble:
	mpiexec -n 120 python3 test_scripts/moist_bubble_3d.py \
		--ncolumns_1d 80 --nlayers 80 --dirname moist_3d_bubble

spatial_convergence_reference:
	python3 test_scripts/linear_sw_gw.py \
		--ncells_per_edge 36 --dt 15.0 --tmax 3600.0 --dumpfreq 240 --stepper trapezium_rule \
		--element_family RTCF --element_order 2; \

spatial_convergence:
	# Triangular elements
	for pair in BDM:1 BDM:2; do \
		family=$${pair%%:*}; \
		order=$${pair##*:}; \
		for ncells in 11 12 13 14 15 16 17 18 19 20; do \
			python3 test_scripts/linear_sw_gw.py \
				--ncells_per_edge $$ncells --dt 60.0 --tmax 3600.0 --dumpfreq 60 --stepper trapezium_rule \
				--element_family $$family --element_order $$order --reference_error \
				--reference_dirname '/data/users/thomas.bendall/gusto_paper_results/results/linear_sw_gw_trapezium_rule_RTCF2_n36_dt15'; \
		done; \
	done

	# Quadrilateral elements
	for pair in RTCF:0 RTCF:1 RTCF:2; do \
		family=$${pair%%:*}; \
		order=$${pair##*:}; \
		for ncells in 15 16 17 18 19 20 21 22 23; do \
			python3 test_scripts/linear_sw_gw.py \
				--ncells_per_edge $$ncells --dt 60.0 --tmax 3600.0 --dumpfreq 60 --stepper trapezium_rule \
				--element_family $$family --element_order $$order --reference_error \
				--reference_dirname '/data/users/thomas.bendall/gusto_paper_results/results/linear_sw_gw_trapezium_rule_RTCF2_n36_dt15'; \
		done; \
	done

linear_convergence:
	# Triangular elements
	for pair in BDM:1 BDM:2; do \
		family=$${pair%%:*}; \
		order=$${pair##*:}; \
		for ncells in 11 12 13 14 15 16 17 18 19 20; do \
			python3 test_scripts/linear_sw_gw.py \
				--ncells_per_edge $$ncells --dt 60.0 --tmax 3600.0 --dumpfreq 60 --stepper trapezium_rule \
				--element_family $$family --element_order $$order; \
		done; \
	done

	# Quadrilateral elements
	for pair in RTCF:0 RTCF:1 RTCF:2; do \
		family=$${pair%%:*}; \
		order=$${pair##*:}; \
		for ncells in 15 16 17 18 19 20 21 22 23; do \
			python3 test_scripts/linear_sw_gw.py \
				--ncells_per_edge $$ncells --dt 60.0 --tmax 3600.0 --dumpfreq 60 --stepper trapezium_rule \
				--element_family $$family --element_order $$order \
		done; \
	done

temporal_convergence:
	for stepper in trapezium_rule ssprk3 rk4 siqn imex_euler leapfrog bdf2; do \
		for pair in 30:360 60:180 120:90; do \
			dt=$${pair%%:*}; \
			dumpfreq=$${pair##*:}; \
			python3 test_scripts/linear_williamson_2.py \
				--ncells_per_edge 16 --dt $$dt --dumpfreq $$dumpfreq \
				--element_family BDM --element_order 1 --stepper $$stepper; \
		done; \
	done
