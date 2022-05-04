all:

tmptest/:
	@mkdir -p tmptest/

tmptest/lowtest.msh: tmptest/
	-@./domain.py -hmesh 200 -o tmptest/lowtest.geo > /dev/null
	-@gmsh -2 tmptest/lowtest.geo > /dev/null

tmptest/slabtest.msh: tmptest/
	-@./domain.py -hmesh 150 -bs 0.0 -o tmptest/slabtest.geo > /dev/null
	-@gmsh -2 tmptest/slabtest.geo > /dev/null

tmptest/margintest.msh: tmptest/
	-@project11/profile.py -hmesh 300 -o tmptest/margintest.geo > /dev/null
	-@gmsh -2 tmptest/margintest.geo > /dev/null

# testing
GMSH = gmsh

rungmshversion_1:
	-@echo "Gmsh version (tested with 4.9.0 and 4.9.3):"
	-@${GMSH} "-version"   # just show gmsh version

rundomain_1: tmptest/
	-@./testit.sh domain.py "-testspew -o tmptest/medtest.geo" 1 1 "generate basic .geo"

runflow_1: tmptest/lowtest.msh
	-@./testit.sh flow.py "-mesh tmptest/lowtest.msh" 1 1 "basic u,p solve with bed step"

runflow_2: tmptest/slabtest.msh
	-@./testit.sh flow.py "-mesh tmptest/slabtest.msh -package Direct" 1 2 "slab u,p solve w direct solver"

runflow_3: tmptest/lowtest.msh
	-@./testit.sh flow.py "-mesh tmptest/lowtest.msh -deltat 20.0 -m 2 -o tmptest/test3.pvd" 1 3 "evolving surface u,p,s solve"

runflow_4: tmptest/lowtest.msh
	-@./testit.sh flow.py "-mesh tmptest/lowtest.msh -save_rank -o tmptest/test4.pvd -s_snes_converged_reason" 2 4 "parallel u,p solve"

runflow_5: tmptest/lowtest.msh
	-@./testit.sh flow.py "-mesh tmptest/lowtest.msh -refine 1 -s_snes_rtol 1.0e-2 -s_ksp_converged_reason -package SchurGMGMassFull -mass_more_reg 1.0" 1 5 "Schur(full)+GMG+nuMass: refine one level"

runflow_6: tmptest/lowtest.msh
	-@./testit.sh flow.py "-mesh tmptest/lowtest.msh -refine 1 -osurface tmptest/test6.png" 1 6 "refine one level and write surface .png"

runflow_7: tmptest/lowtest.msh
	-@./testit.sh flow.py "-mesh tmptest/lowtest.msh -sequence 1 -s_snes_converged_reason -package SchurGMGMass" 1 7 "Schur+GMG+nuMass: grid-sequence one level"

runflow_8: tmptest/margintest.msh
	-@./testit.sh flow.py "-mesh tmptest/margintest.msh -alpha 0.0" 1 8 "basic solve for project11 margin shape"

test: rungmshversion_1 rundomain_1 runflow_1 runflow_2 runflow_3 runflow_4 runflow_5 runflow_6 runflow_7 runflow_8

.PHONY: clean rungmshversion_1 rundomain_1 runflow_1 runflow_2 runflow_3 runflow_4 runflow_5 runflow_6 runflow_7 runflow_8 test

clean:    # aggressive clean in current directory and doc/,study/,project11/
	rm -f *~ *.pyc *.pvd *.vtu *.pvtu *.geo *.msh *.png
	rm -rf __pycache__/ tmptest/
	(cd doc/; ${MAKE} clean)
	(cd study/; ${MAKE} clean)
	(cd project11/; ${MAKE} clean)
