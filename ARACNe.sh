sample=sample_name
cd /path/ARACNe
#mkdir ARACNe
mkdir "$sample"
mkdir "$sample"/bootstraps_tfs
mkdir "$sample"/bootstraps_cotfs
mkdir "$sample"/bootstraps_sig
mkdir "$sample"/bootstraps_surface

cd /path/ARACNe-AP/dist
#tfs
#calculate threshold
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -e /path/ARACNe/"$sample"_cpm.tsv \
   -o /path/ARACNe/"$sample"/bootstraps_tfs \
   --tfs /path/TF/tfs-ensembl.txt \
   --pvalue 1E-8 --seed 1 --calculateThreshold;
#run bootstraps
for i in {1..100}
do
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -e /path/ARACNe/"$sample"_cpm.tsv \
   -o /path/ARACNe/"$sample"/bootstraps_tfs \
   --tfs /path/TF/tfs-ensembl.txt \
   --pvalue 1E-8 --seed $i
done;
#consolidate
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -o /path/ARACNe/"$sample"/bootstraps_tfs \
   --consolidate;

#cotfs
#calculate threshold
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -e /path/ARACNe/"$sample"_cpm.tsv \
   -o /path/ARACNe/"$sample"/bootstraps_cotfs \
   --tfs /path/TF/cotfs-ensembl.txt \
   --pvalue 1E-8 --seed 1 --calculateThreshold;
#run bootstraps
for i in {1..100}
do
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -e /path/ARACNe/"$sample"_cpm.tsv \
   -o /path/ARACNe/"$sample"/bootstraps_cotfs \
   --tfs /path/TF/cotfs-ensembl.txt \
   --pvalue 1E-8 --seed $i
done;
#consolidate
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -o /path/ARACNe/"$sample"/bootstraps_cotfs \
   --consolidate;

 #sig
#calculate threshold
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -e /path/ARACNe/"$sample"_cpm.tsv \
   -o /path/ARACNe/"$sample"/bootstraps_sig \
   --tfs /path/TF/sig-ensembl.txt \
   --pvalue 1E-8 --seed 1 --calculateThreshold;
#run bootstraps
for i in {1..100}
do
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -e /path/ARACNe/"$sample"_cpm.tsv \
   -o /path/ARACNe/"$sample"/bootstraps_sig \
   --tfs /path/TF/sig-ensembl.txt \
   --pvalue 1E-8 --seed $i
done;
#consolidate
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -o /path/ARACNe/"$sample"/bootstraps_sig \
   --consolidate;

 #surface
#calculate threshold
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -e /path/ARACNe/"$sample"_cpm.tsv \
   -o /path/ARACNe/"$sample"/bootstraps_surface \
   --tfs /path/TF/surface-ensembl.txt \
   --pvalue 1E-8 --seed 1 --calculateThreshold;
#run bootstraps
for i in {1..100}
do
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -e /path/ARACNe/"$sample"_cpm.tsv \
   -o /path/ARACNe/"$sample"/bootstraps_surface \
   --tfs /path/TF/surface-ensembl.txt \
   --pvalue 1E-8 --seed $i
done;
#consolidate
~/anaconda3/bin/java -Xmx5G -jar aracne.jar \
   -o /path/ARACNe/"$sample"/bootstraps_surface \
   --consolidate;

mkdir network

cp bootstraps_tfs/network.txt network/network_tfs.txt
cp bootstraps_cotfs/network.txt network/network_cotfs.txt
cp bootstraps_sig/network.txt network/network_sig.txt
cp bootstraps_surface/network.txt network/network_surface.txt

cd network
awk 'NR>1' network_tfs.txt > tmp && mv tmp network_tfs.txt
awk 'NR>1' network_cotfs.txt > tmp && mv tmp network_cotfs.txt
awk 'NR>1' network_sig.txt > tmp && mv tmp network_sig.txt
awk 'NR>1' network_surface.txt > tmp && mv tmp network_surface.txt
cat network_tfs.txt network_cotfs.txt network_sig.txt network_surface.txt >> network_final.txt

