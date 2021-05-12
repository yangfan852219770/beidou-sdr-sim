#Makefile
objs = lime_beidou.o beidou_sim.o util.o
test: $(objs)
	gcc -o $@ $^ -lm -lpthread -lLimeSuite
beidou_sim.o: beidou_sim.h lime_beidou.h util.h
lime_beidou.o: lime_beidou.h
util.o: util.h
.PHONY : clean
clean : 
	-rm *.o
