#---INTERFACE TO CONTROLLER
#-------------------------------------------------------------------------------------------------------------

#---valid function names from the python script
TARGETS := $(shell perl -n -e '@parts = /^def\s+[a-z,_]+/g; $$\ = "\n"; print for @parts;' \
	controller.py | awk '{print $$2}')

#---collect arguments
RUN_ARGS := $(wordlist 1,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
$(eval $(RUN_ARGS):;@:)

#---hack to always run makefile interfaced with python
scripts=controller.py
$(shell touch $(scripts))
checkfile=.pipeline_up_to_date

#---targets
$(checkfile): $(scripts)
	touch $(checkfile)
	python controller.py ${RUN_ARGS} ${MAKEFLAGS}

#---default and arbitrary make targets
default: $(checkfile)
$(TARGETS): $(checkfile)

#---git push
push:
	bash amx/script-gitpush.sh ${RUN_ARGS}

#---quickly view simulations in VMD
view:
	touch $(checkfile)
	python amx/script-vmd-viewer.py ${RUN_ARGS} ${MAKEFLAGS}
