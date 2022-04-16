# Makefile for MATH 214 Capstone Project

CXX ?= g++
CXXFLAGS ?= -Wall -Werror -pedantic --std=c++17 -g


spectral_clustering.exe: spectral_clustering.cpp
	$(CXX) $(CXXFLAGS) spectral_clustering.cpp -o spectral_clustering.exe