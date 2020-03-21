build-run:
	make rebuild
	make run

rebuild:
	cmake -B build -DCMAKE_BUILD_TYPE=Release 
	make -C build -j

run:
	./build/Practical3_bin ./data/ bunny-scene.txt emptyconstraints.txt
