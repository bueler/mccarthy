
lowtest.msh:
	-@./gendomain.py -hmesh 200 -o lowtest.geo > /dev/null
	-@gmsh -2 lowtest.geo > /dev/null

slabtest.msh:
	-@./gendomain.py -hmesh 150 -bs 0.0 -o slabtest.geo > /dev/null
	-@gmsh -2 slabtest.geo > /dev/null

# testing
rungendomain_1:
	-@./testit.sh gendomain.py "-testspew -o medtest.geo" 1 1 "generate basic .geo"

runflow_1: lowtest.msh
	-@./testit.sh flow.py "lowtest.msh" 1 1 "basic u,p solve with bed step"

runflow_2: slabtest.msh
	-@./testit.sh flow.py "slabtest.msh" 1 2 "slab bed u,p solve"

runflow_3: lowtest.msh
	-@./testit.sh flow.py "-deltat 20.0 -m 2 -o test3.pvd lowtest.msh" 1 3 "evolving surface u,p,h solve"

runflow_4: lowtest.msh
	-@./testit.sh flow.py "-o test4.pvd lowtest.msh -s_snes_converged_reason" 2 4 "parallel u,p solve"

test: rungendomain_1 runflow_1 runflow_2 runflow_3 runflow_4

.PHONY: clean rungendomain_1 runflow_1 runflow_2 runflow_3 runflow_4 test

clean:
	@rm -f *~ *.pyc *.pvd *.vtu *.pvtu *.geo *.msh
	@rm -rf __pycache__/

