all:parallel1

serial:main.cpp
	g++ main.cpp
	./a.out


parallel:main.cpp
	./c.sh 4

parallel1:main.cpp
	mpic++ -DUSE_MPI main.cpp
	mpirun -n 4 ./a.out

clean:
#	rm -f ./*.txt
	rm -f ./*.figdata
	rm -f *.out

rmdata:
	rm -f ./Data*.txt
	rm -f ./Run-conf*.txt
	rm -f *.out
	rm -f ID.txt

	
	
