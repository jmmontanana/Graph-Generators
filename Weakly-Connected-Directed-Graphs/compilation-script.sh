#!/bin/bash

gcc -o duo duo.c write_r.c grafos.c tamano_com.c -lm -Wall;
gcc -o dbo dbo.c write_r.c grafos.c tamano_com.c -lm -Wall;

