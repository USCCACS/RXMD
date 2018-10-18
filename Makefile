default:
	@echo " rxmd package "
	@echo "  make <OPTIONS>"
	@echo
	@echo " OPTIONS:"
	@echo "   all"
	@echo "   init"
	@echo "   rxmd"
	@echo "   clean"

.PHONY: clean all rxmd init
# both init and rxmd
all:
	cd init && $(MAKE)
	cd src && $(MAKE)
rxmd:
	cd src && $(MAKE)
init:
	cd init && $(MAKE)
clean:
	cd init && $(MAKE) clean
	cd src && $(MAKE) clean
