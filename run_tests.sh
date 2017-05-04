#!/bin/bash
# Script to run tests for paper in JASSS

if [[ "$#" -eq 0 ]]
then
    echo "Running all tests."
    TIME_ALGS=1
    TIME_CSPM_AGENTS=1
    TIME_K=1
    ES5K=1
    EA5K=1
    ES20K=1
    EK=1
    EC=1
fi

while (( "$#" ));
do
    case $1 in
        -t|--timealgs)
            TIME_ALGS=1
            ;;
        -a|--timeagents)
            TIME_CSPM_AGENTS=1
            ;;
        -k|--timek)
            TIME_K=1
            ;;
        -es5k)
            ES5K=1
            ;;
        -ea5k)
            EA5K=1
            ;;
        -es20k)
            ES20K=1
            ;;
        -ek)
            EK=1
            ;;
        -ec)
            EC=1
            ;;
        *)
            echo "Unknown argument: " $1
            echo "See README for usage"
            exit
            # unknown option
            ;;
    esac
    shift # past argument or value
done

# TIMING TESTS

if [[ "$TIME_ALGS" -eq 1 ]]
then
    echo "Running timing tests for algorithms"
    echo "*******************" >> tables.tex
    echo "% Timing: All algorithms" >tables.tex
    make release
    ./partnermatch -t -h -n 20000 -k 200 -c 100 -r 10 -f 0 -i 20 -s $RANDOM >timing_tests.csv
    Rscript analyzeAlgorithmTimes.R >>tables.tex
fi

if [[ "$TIME_CSPM_AGENTS" -eq 1 ]]
then
    echo "Running timing tests for CSPM on increasing agents"
    echo "% Timing: Increasing agents" >>tables.tex
    make release
    ./partnermatch -t -h -d 50000 -n 50000 -k 200 -c 100 -r 10 -d 50000 -i 1 -s -a C $RANDOM >timing_agent_increase_tests.csv
    ./partnermatch -t -n 100000 -k 200 -c 100 -r 10 -d 100000 -i 1 -a C -s $RANDOM >>timing_agent_increase_tests.csv
    ./partnermatch -t -n 500000 -k 200 -c 100 -r 10 -d 500000 -i 1 -a C -s $RANDOM >>timing_agent_increase_tests.csv
    ./partnermatch -t -n 1000000 -k 200 -c 100 -r 10 -d 1000000 -i 1 -a C -s $RANDOM >>timing_agent_increase_tests.csv
    ./partnermatch -t -n 5000000 -k 200 -c 100 -r 10 -d 5000000 -i 1 -a C -s $RANDOM >>timing_agent_increase_tests.csv
    Rscript analyzeIncreasingAgentCSPMTimes.R >>tables.tex
fi

if [[ "$TIME_K" -eq 1 ]]
then
    echo "Running timing tests for CSPM on k"
    echo "*******************" >> tables.tex
    echo "% Timing: CSPM on k" >>tables.tex
    make release
    ./partnermatch -t -h -n 500000 -k 50 -k 100 -r 10 -i 1 -a C -s $RANDOM >timing_k_tests.csv
    ./partnermatch -t  -n 500000 -k 100 -c 100 -r 10 -i 1 -a C -s $RANDOM >>timing_k_tests.csv
    ./partnermatch -t  -n 500000 -k 150 -c 100 -r 10 -i 1 -a C -s $RANDOM >>timing_k_tests.csv
    ./partnermatch -t  -n 500000 -k 200 -c 100 -r 10 -i 1 -a C -s $RANDOM >>timing_k_tests.csv
    ./partnermatch -t  -n 500000 -k 250 -c 100 -r 10 -i 1 -a C -s $RANDOM >>timing_k_tests.csv
    ./partnermatch -t  -n 500000 -k 300 -c 100 -r 10 -i 1 -a C -s $RANDOM >>timing_k_tests.csv
    ./partnermatch -t  -n 500000 -k 350 -c 100 -r 10 -i 1 -a C -s $RANDOM >>timing_k_tests.csv
    ./partnermatch -t  -n 500000 -k 400 -c 100 -r 10 -i 1 -a C -s $RANDOM >>timing_k_tests.csv
    ./partnermatch -t  -n 500000 -k 450 -c 100 -r 10 -i 1 -a C -s $RANDOM >>timing_k_tests.csv
    ./partnermatch -t  -n 500000 -k 500 -c 100 -r 10 -i 1 -a C -s $RANDOM >>timing_k_tests.csv
    Rscript analyzeKSpeeds.R >>tables.tex
fi

# END OF TIMING TESTS

# BEGINNING OF EFFECTIVENESS TESTS

if [[ "$ES5K" -eq 1 ]]
then
    echo "Running STI effectiveness tests:5000"
    echo "*******************" >> tables.tex
    echo "% Effectiveness STI 5000" >>tables.tex
    rm output5kS_??.csv
    make release
    ./partnermatch -h -d 1 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_01.csv &
    ./partnermatch -d 2 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_02.csv &
    ./partnermatch -d 3 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_03.csv &
    ./partnermatch -d 4 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_04.csv &
    ./partnermatch -d 5 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_05.csv &
    ./partnermatch -d 6 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_06.csv &
    ./partnermatch -d 7 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_07.csv &
    ./partnermatch -d 8 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_08.csv &
    ./partnermatch -d 9 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_09.csv &
    ./partnermatch -d 10 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_10.csv &
    ./partnermatch -d 11 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_11.csv &
    ./partnermatch -d 12 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_12.csv &
    ./partnermatch -d 13 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_13.csv &
    ./partnermatch -d 14 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_14.csv &
    ./partnermatch -d 15 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_15.csv &
    ./partnermatch -d 16 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_16.csv &
    ./partnermatch -d 17 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_17.csv &
    ./partnermatch -d 18 -n 5000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM -b >output5kS_18.csv &
    wait
    echo "Analyzing STI effectiveness tests:5000"
    rm output5kS.csv
    cat output5kS_??.csv >output5kS.csv
    cp output5kS.csv tmp_analysis.csv
    Rscript analyzeEffectiveness.R >>tables.tex
fi


if [[ "$EA5K" -eq 1 ]]
then
   echo "Running ATTRACTOR effectiveness tests:5000"
   echo "*******************" >> tables.tex
   echo "% Effectiveness ATTRACTOR 5000" >>tables.tex
   make release_attract
   ./partnermatch -h -d 0.00 -n 5000 -r 20 -i 1 -k 200 -c 100 -A 0 -R 1 -b > output5k_ra_000.csv &
   ./partnermatch -h -d 0.25 -n 5000 -r 20 -i 1 -k 200 -c 100 -A 0.25 -R 0.75 -b > output5k_ra_025.csv &
   ./partnermatch -h -d 0.50 -n 5000 -r 20 -i 1 -k 200 -c 100 -A 0.50 -R 0.50 -b > output5k_ra_050.csv &
   ./partnermatch -h -d 0.75 -n 5000 -r 20 -i 1 -k 200 -c 100 -A 0.75 -R 0.25 -b > output5k_ra_075.csv &
   ./partnermatch -h -d 1.00 -n 5000 -r 20 -i 1 -k 200 -c 100 -A 1 -R 0 -b > output5k_ra_100.csv &
   wait
   echo "****************" >>tables.tex
   echo "Attractor = 0.00" >>tables.tex
   cp output5k_ra_000.csv tmp_analysis.csv
   Rscript analyzeEffectiveness.R >>tables.tex

   echo "****************" >>tables.tex
   echo "Attractor = 0.25" >>tables.tex
   cp output5k_ra_025.csv tmp_analysis.csv
   Rscript analyzeEffectiveness.R >>tables.tex

   echo "****************" >>tables.tex
   echo "Attractor = 0.50" >>tables.tex
   cp output5k_ra_050.csv tmp_analysis.csv
   Rscript analyzeEffectiveness.R >>tables.tex

   echo "****************" >>tables.tex
   echo "Attractor = 0.75" >>tables.tex
   cp output5k_ra_075.csv tmp_analysis.csv
   Rscript analyzeEffectiveness.R >>tables.tex

   echo "****************" >>tables.tex
   echo "Attractor = 1.00" >>tables.tex
   cp output5k_ra_100.csv tmp_analysis.csv
   Rscript analyzeEffectiveness.R >>tables.tex

   echo "****************" >>tables.tex
   echo "Attractor tables" >>tables.tex
   Rscript analyzeAttractor.R >>tables.tex
fi


if [[ "$ES20K" -eq 1 ]]
then
    echo "Running STI effectiveness tests:20000"
    echo "*******************" >> tables.tex
    echo "% Effectiveness STI 20000" >>tables.tex
    rm output20kS_??.csv
    make release
    ./partnermatch -h -d 1 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_01.csv &
    ./partnermatch -d 2 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_02.csv &
    ./partnermatch -d 3 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_03.csv &
    ./partnermatch -d 4 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_04.csv &
    ./partnermatch -d 5 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_05.csv &
    ./partnermatch -d 6 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_06.csv &
    ./partnermatch -d 7 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_07.csv &
    ./partnermatch -d 8 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_08.csv &
    ./partnermatch -d 9 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_09.csv &
    ./partnermatch -d 10 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_10.csv &
    ./partnermatch -d 11 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_11.csv &
    ./partnermatch -d 12 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_12.csv &
    ./partnermatch -d 13 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_13.csv &
    ./partnermatch -d 14 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_14.csv &
    ./partnermatch -d 15 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_15.csv &
    ./partnermatch -d 16 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_16.csv &
    ./partnermatch -d 17 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_17.csv &
    ./partnermatch -d 18 -n 20000 -r 2 -i 20 -k 200 -c 100 -s $RANDOM  >output20kS_18.csv &
    wait
    echo "Analyzing STI effectiveness tests:20000"
    rm output20kS.csv
    cat output20kS_??.csv >output20kS.csv
    cp output20kS.csv tmp_analysis.csv
    Rscript analyzeEffectiveness.R >>tables.tex
    echo "************" >> tables.tex
    echo "Best/worst for 20k" >> tables.tex
    Rscript bestWorst.R >>tables.tex
fi


if [[ "$EK" -eq 1 ]]
then
    echo "Running K tests"
    rm output20k_k????.csv
    make release
    ./partnermatch -h -n 20000 -a C -r 40 -i 1 -k 50 -c 100 -s $RANDOM >output20k_k0050.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 100 -c 100 -s $RANDOM >output20k_k0100.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 150 -c 100 -s $RANDOM >output20k_k0150.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 200 -c 100 -s $RANDOM >output20k_k0200.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 250 -c 100 -s $RANDOM >output20k_k0250.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 300 -c 100 -s $RANDOM >output20k_k0300.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 350 -c 100 -s $RANDOM >output20k_k0350.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 400 -c 100 -s $RANDOM >output20k_k0400.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 450 -c 100 -s $RANDOM >output20k_k0450.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 500 -c 100 -s $RANDOM >output20k_k0500.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 550 -c 100 -s $RANDOM >output20k_k0550.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 600 -c 100 -s $RANDOM >output20k_k0600.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 650 -c 100 -s $RANDOM >output20k_k0650.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 700 -c 100 -s $RANDOM >output20k_k0700.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 750 -c 100 -s $RANDOM >output20k_k0750.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 800 -c 100 -s $RANDOM >output20k_k0800.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 850 -c 100 -s $RANDOM >output20k_k0850.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 900 -c 100 -s $RANDOM >output20k_k0900.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 950 -c 100 -s $RANDOM >output20k_k0950.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -k 1000 -c 100 -s $RANDOM >output20k_k1000.csv &
    wait
    cat output20k_k????.csv >output20k_k.csv
    Rscript analyzeK.R
fi

if [[ "$EC" -eq 1 ]]
then
    echo "Running C tests"
    rm output20k_c????.csv
    make release
    ./partnermatch -h -n 20000 -a C -r 40 -i 1 -c 1 -k 200 -s $RANDOM >output20k_c0001.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 50 -k 200 -s $RANDOM >output20k_c0050.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 100 -k 200 -s $RANDOM >output20k_c0100.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 150 -k 200 -s $RANDOM >output20k_c0150.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 200 -k 200 -s $RANDOM >output20k_c0200.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 250 -k 200 -s $RANDOM >output20k_c0250.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 300 -k 200 -s $RANDOM >output20k_c0300.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 350 -k 200 -s $RANDOM >output20k_c0350.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 400 -k 200 -s $RANDOM >output20k_c0400.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 450 -k 200 -s $RANDOM >output20k_c0450.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 500 -k 200 -s $RANDOM >output20k_c0500.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 550 -k 200 -s $RANDOM >output20k_c0550.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 600 -k 200 -s $RANDOM >output20k_c0600.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 650 -k 200 -s $RANDOM >output20k_c0650.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 700 -k 200 -s $RANDOM >output20k_c0700.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 750 -k 200 -s $RANDOM >output20k_c0750.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 800 -k 200 -s $RANDOM >output20k_c0800.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 850 -k 200 -s $RANDOM >output20k_c0850.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 900 -k 200 -s $RANDOM >output20k_c0900.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 950 -k 200 -s $RANDOM >output20k_c0950.csv &
    ./partnermatch -n 20000 -a C -r 40 -i 1 -c 1000 -k 200 -s $RANDOM >output20k_c1000.csv &
    wait
    cat output20k_c????.csv >output20k_c.csv
    Rscript analyzeC.R
fi
