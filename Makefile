### make changes accordingly ###
CC       = gcc
CPP      = g++
CLINKER  = gcc
CCLINKER = g++
MAKE     = make --no-print-directory
SHELL    = /bin/sh
CFLAGS		= -Wall -pedantic 
OPTI            = -O3
#OPTI = -pg -g # for profiling
LDFLAGS	= 
INCLUDES	= 
LIBS	= -lm -lgsl -lgslcblas
OBJS	= GlobalModel.o CityModel.o
DEFINES = -DVERBOSE 

default: global

global: $(OBJS) Makefile driver.o
	$(CCLINKER) -o global driver.o $(OBJS) $(LDFLAGS) $(LIBS)

mixedglobal: $(OBJS) Makefile mixingdriver.o
	$(CCLINKER) -o mixedglobal mixingdriver.o $(OBJS) $(LDFLAGS) $(LIBS)

yearlyglobal: $(OBJS) Makefile yearlydriver.o
	$(CCLINKER) -o yearlyglobal yearlydriver.o $(OBJS) $(LDFLAGS) $(LIBS)

CityModel.o: CityModel.cpp CityModel.h GlobalModel.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c CityModel.cpp

GlobalModel.o: GlobalModel.cpp CityModel.h GlobalModel.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c GlobalModel.cpp

driver.o: driver.cpp CityModel.h GlobalModel.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c driver.cpp

yearlydriver.o: yearlydriver.cpp CityModel.h GlobalModel.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c yearlydriver.cpp

%.o: %.cpp CityModel.h GlobalModel.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c $<

zip: *.cpp *.h Makefile README gpl.txt DataFiles/*
	zip resistancesimulation.zip README gpl.txt Makefile *.cpp *.h DataFiles/population_321_age.txt DataFiles/travel_321.txt DataFiles/seasonality_321.csv

emacs:
	emacs Makefile *.h *.cpp *R ../manuscript/resistance.tex &

clean:
	rm -f *.o global yearlyglobal mixedglobal *~
