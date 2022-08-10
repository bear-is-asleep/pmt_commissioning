for i in {1..5}
do
  python3.9 make_pkls.py $i gun0
  #python3.9 make_pkls.py $i ew
  #python3.9 make_pkls.py $i fb
  echo $i
done