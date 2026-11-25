#!/bin/sh
# Step 13b gate-length study (2026-11-18), Step 13c sweep extension (2026-11-23), Step 13d snapshots (2026-11-24).
# Run inside the Device 2 folder (needs sdevice_n5.par and, for ext, the Step 13b meshes NAME_msh.tdr).
#   sh run_lg.sh mesh               build all 12 meshes  (sde, about a minute each)
#   sh run_lg.sh run g050           run one device       (sdevice, in the background)
#   sh run_lg.sh queue g100 ...     run several, one after another, in the background
#   sh run_lg.sh ext g500           Step 13c: +0.8 -> +2 V extension of one device (lg_ext_des.cmd -> NAMEx_*)
#   sh run_lg.sh extqueue g500 ...  Step 13c: several extensions, one after another, in the background
#   sh run_lg.sh snap a200          Step 13d: extension with band/density snapshots (lg_snap_des.cmd -> NAMEs_*)
# Names: g = graded [P], a = abrupt control [ASSUME]; 050 = L_G 0.05 um ... 700 = 0.70 um.
ALL="g050 g100 g200 g350 g500 g700 a050 a100 a200 a350 a500 a700"
lg()  { echo "0.$(echo "$1" | cut -c2-4)" | sed 's/^0\.\(.\)\(..\)$/0.\1\2/' ; }
gr()  { case "$1" in g*) echo 1 ;; a*) echo 0 ;; esac ; }
mk()  { # $1 name -> ${1}_dvs.cmd, ${1}_full.cmd
  L=$(lg "$1"); G=$(gr "$1")
  sed -e "s/@LG@/$L/g" -e "s/@GRADED@/$G/g" -e "s/n@node@/$1/g" sde_lg.cmd > "$1"_dvs.cmd
  sed -e "s/@NAME@/$1/g" -e "s/@LG@/$L/g" -e "s/@GRADED@/$G/g" lg_des.cmd > "$1"_full.cmd ; }
mkx() { # $1 name -> ${1}x_full.cmd (needs ${1}_msh.tdr from Step 13b)
  [ -f "$1"_msh.tdr ] || { echo "missing $1_msh.tdr (run: sh run_lg.sh mesh)"; return 1; }
  sed -e "s/@NAME@/$1/g" lg_ext_des.cmd > "$1"x_full.cmd ; }
mks() { # $1 name -> ${1}s_full.cmd (needs ${1}_msh.tdr from Step 13b)
  [ -f "$1"_msh.tdr ] || { echo "missing $1_msh.tdr (run: sh run_lg.sh mesh)"; return 1; }
  sed -e "s/@NAME@/$1/g" lg_snap_des.cmd > "$1"s_full.cmd ; }
case "$1" in
  mesh)  for n in $ALL; do mk $n; sde -e -l ${n}_dvs.cmd > ${n}_dvs.out 2>&1
           echo "$n: L_G $(lg $n) um, graded $(gr $n): $(grep -c 'Meshing successful' ${n}_dvs.out) x 'Meshing successful'"; done ;;
  run)   mk $2; nohup sdevice ${2}_full.cmd > ${2}_des.out 2>&1 &
         echo "started $2 (pid $!)" ;;
  queue) shift; for n in "$@"; do mk $n; done
         nohup sh -c 'for n in "$@"; do sdevice ${n}_full.cmd > ${n}_des.out 2>&1; done' sh "$@" > lg_queue.out 2>&1 &
         echo "queued $* (pid $!)" ;;
  ext)   mkx $2 || exit 1; nohup sdevice ${2}x_full.cmd > ${2}x_des.out 2>&1 &
         echo "started $2 extension (pid $!)" ;;
  extqueue) shift; for n in "$@"; do mkx $n || exit 1; done
         nohup sh -c 'for n in "$@"; do sdevice ${n}x_full.cmd > ${n}x_des.out 2>&1; done' sh "$@" > lg_extqueue_$1.out 2>&1 &
         echo "queued extensions $* (pid $!)" ;;
  snap)  mks $2 || exit 1; nohup sdevice ${2}s_full.cmd > ${2}s_des.out 2>&1 &
         echo "started $2 snapshot run (pid $!)" ;;
  *)     echo "usage: sh run_lg.sh mesh | run NAME | queue NAME... | ext NAME | extqueue NAME... | snap NAME" ;;
esac
