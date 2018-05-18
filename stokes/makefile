
lowtest.msh:
	-@./gendomain.py -hmesh 200 -o lowtest.geo > /dev/null
	-@gmsh -2 lowtest.geo > /dev/null

slabtest.msh:
	-@./gendomain.py -hmesh 150 -bs 0.0 -o slabtest.geo > /dev/null
	-@gmsh -2 slabtest.geo > /dev/null

# testing
rungendomain_1:
	-@./testit.sh gendomain.py "-testspew -o medtest.geo" 1 1

runflow_1: lowtest.msh
	-@./testit.sh flow.py "lowtest.msh" 1 1

runflow_2: slabtest.msh
	-@./testit.sh flow.py "slabtest.msh" 1 2

test: rungendomain_1 runflow_1 runflow_2

.PHONY: clean rungendomain_1 runflow_1 runflow_2 test

clean:
	@rm -f *~ *.pyc *.pvd *.vtu *.pvtu *.geo *.msh
	@rm -rf __pycache__/

