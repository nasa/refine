#!/usr/bin/env bash

bamg -g square_g.msh -o square_0.msh -thetaquad 75

ref translate square_0.msh square_0.tec
ref translate square_0.msh square_0.b8.ugrid

