all:

tmptest/:
	mkdir -p tmptest/

tmptest/lowtest.msh: tmptest/
	-@./gendomain.py -hmesh 200 -o tmptest/lowtest.geo > /dev/null
	-@gmsh -2 tmptest/lowtest.geo > /dev/null

tmptest/slabtest.msh: tmptest/
	-@./gendomain.py -hmesh 150 -bs 0.0 -o tmptest/slabtest.geo > /dev/null
	-@gmsh -2 tmptest/slabtest.geo > /dev/null

# testing
rungendomain_1:
	-@./testit.sh gendomain.py "-testspew -o tmptest/medtest.geo" 1 1 "generate basic .geo"

runflow_1: tmptest/lowtest.msh
	-@./testit.sh flow.py "tmptest/lowtest.msh" 1 1 "basic u,p solve with bed step"

runflow_2: tmptest/slabtest.msh
	-@./testit.sh flow.py "tmptest/slabtest.msh" 1 2 "slab u,p solve"

runflow_3: tmptest/lowtest.msh
	-@./testit.sh flow.py "-deltat 20.0 -m 2 -o tmptest/test3.pvd tmptest/lowtest.msh" 1 3 "evolving surface u,p,h solve"

runflow_4: tmptest/lowtest.msh
	-@./testit.sh flow.py "-o tmptest/test4.pvd tmptest/lowtest.msh -s_snes_converged_reason" 2 4 "parallel u,p solve"

test: rungendomain_1 runflow_1 runflow_2 runflow_3 runflow_4

.PHONY: clean rungendomain_1 runflow_1 runflow_2 runflow_3 runflow_4 test

clean:
	rm -f *~ *.pyc *.pvd *.vtu *.pvtu *.geo *.msh
	rm -rf __pycache__/ tmptest/
	(cd doc/; ${MAKE} clean)

