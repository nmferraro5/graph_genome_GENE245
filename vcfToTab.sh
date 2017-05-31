gzcat ALL.chr2.noCom.vcf.gz | cut -f1-5,8 | bgzip -c > ALL.chr2.afs.vcf.gz
gzcat ALL.chr2.afs.vcf.gz | head -10 | sed 's/;/       /g'
Remove descr= pattern
cat ALL.chr2.afs2.txt | sed 's/..=//g' | sed 's/..._//g' > ALL.chr2.afs3.txt
Add column names via vim
