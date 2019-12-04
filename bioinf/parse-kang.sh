# shell script to rearrange output from G-quad finding code into description of reverting vs non-reverting pairs in Kang et al. data

# label actual number of Cs in a particular borderline case (8_H44a) 
sed 's/0 0 0/2 7 4/g' kang-fasta.fasta-out.txt > tmp

# rearrange outputs for parsing in R
awk 'BEGIN{
  printf("label\tl_300_diff\tl_16180_diff\tl_16345_diff\tDid_A_win\tref\n");
  n = 0;
}{
  c1 = $3+$5; c2 = $7+$9; c3 = $11+$13;
  if(n % 2 == 1)
  {
    printf("%i\t%i\t%i\t%i\t%i\t%i\n", n/2+1, oldc1-c1, oldc2-c2, oldc3-c3, (n/2+1 <= 4.5 ? 1 : 0), n/2+1);
  }
  oldc1 = c1; oldc2 = c2; oldc3 = c3;
  n++;
}' tmp > kang-competition.txt
