for i in {1..5}
do 
  #hadd -f ew_$i.root /pnfs/sbnd/scratch/users/brindenc/v09_43_00/ew_filtered/hitdumper_ew/*_$i/hitdumper_tree.root
  #hadd -f fb_$i.root /pnfs/sbnd/scratch/users/brindenc/v09_43_00/fb_filtered/hitdumper_fb/*_$i/hitdumper_tree.root
  hadd -f gun0_$i.root /pnfs/sbnd/scratch/users/brindenc/v09_54_00/gun_gun/gun0_hitdumper/*_$i/hitdumper_tree.root
done
