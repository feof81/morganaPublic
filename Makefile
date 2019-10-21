###### Main Rules
all: 
	cd src && $(MAKE) all && cd ..

tau:
	cd src && $(MAKE) all -f ./makeTau && cd ..
	
clean:
	cd src && $(MAKE) clean && cd ..


