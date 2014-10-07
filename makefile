#---INTERFACE TO CONTROLLER
#-------------------------------------------------------------------------------------------------------------

#---extract arguments for "render" function
ifeq (render,$(firstword $(MAKECMDGOALS)))
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(RUN_ARGS):;@:)
endif

help:
	./controller
aamd-bilayer: 
	./controller make aamd-bilayer	
cgmd-bilayer: 
	./controller make cgmd-bilayer	
aamd-protein: 
	./controller make aamd-protein	
cgmd-protein: 
	./controller make cgmd-protein
clean: 
	./controller reset
reset: 
	./controller reset
render:
	./controller render $(RUN_ARGS)
lastframe:
	./amx/script-lastframe.py

#---since docs is a folder we depend on the code from whence docs is generated
docs: amx
	./controller make_docs
	
#---quickly testing development code
dev:
	echo -e "y\ny\n" | ./controller reset
	./controller make aamd-protein
	./script-master-aamd-protein
