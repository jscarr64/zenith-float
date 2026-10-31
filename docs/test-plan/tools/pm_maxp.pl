# Harper pmrev12 A2 (Paul ruling 2026-10-30 8:06 PM CT): PM_MAX_P resolves from zf's MAX_PRECISION in the source the matrix is pinned to.
# No MAX value is written in any tool: the value comes from git show <commit>:<path>, the expected M32 value from inputs/pm-constants.csv.
use strict; use warnings; use Math::BigInt;
sub pm_maxp { my ($T,$commit,$pcf)=@_;
  my %pc; open my $ch,'<',$pcf or die "$pcf\n"; <$ch>; while(<$ch>){ chomp; next if /^\s*$/; my @c=split /,/,$_,-1; $pc{$c[0]}=[@c] }
  my $e=$pc{PM_MAX_P} or die "pm-constants.csv lacks PM_MAX_P\n";
  my ($path,$const)=$e->[1]=~/^src:([^:]+):(\w+)$/ or die "PM_MAX_P value '$e->[1]' must be src:<path>:<CONST> (never a literal)\n";
  my $txt=`git -C "$T" show $commit:$path`; die "git show $commit:$path failed\n" if $?;
  my @L=split /\n/,$txt,-1; my ($ln)=grep { $L[$_]=~/^pub const \Q$const\E: Precision = / } 0..$#L; die "$const not found in $path at $commit\n" unless defined $ln;
  my ($v,$form);
  if($L[$ln]=~/= ([0-9_]+);\s*$/){ ($v=$1)=~s/_//g; $v=Math::BigInt->new($v); $form='integer literal' }
  elsif($L[$ln]=~/= \{\s*$/ && join("\n",@L[$ln+1..$ln+4])=~/\A\s*let by_exponent: u64 = 1 << (\d+);\n\s*let by_size = \(usize::MAX \/ (\d+)\) as u64;\n\s*let lim = if by_exponent < by_size \{ by_exponent \} else \{ by_size \};\n\s*\(lim - 1\) as Precision;?\z/){
    my $be=Math::BigInt->new(2)->bpow($1); my $bs=Math::BigInt->new(2)->bpow(64)->bsub(1)->bdiv($2); my $lim=$be->bcmp($bs)<0?$be:$bs; $v=$lim->copy->bsub(1); $form="min(2^$1, usize::MAX/$2) - 1, evaluated for 64-bit usize (P6)" }
  else { die "UNRESOLVED: $path:".($ln+1)." $const has a form tools/pm_maxp.pl does not evaluate\n" }
  my $sha=`git -C "$T" show $commit:$path | sha256sum`; $sha=~s/\s.*//s;
  my ($xv,$xc,$xl)=$e->[3]=~/^expected=(\d+)\@([0-9a-f]{40}):(\S+)$/ or die "PM_MAX_P status must be expected=<int>\@<sha40>:<file:line>\n";
  my $m32=system("git","-C",$T,"merge-base","--is-ancestor",$xc,$commit)==0?1:0;
  return {value=>"$v",form=>$form,file=>$path,line=>$ln+1,sha=>$sha,commit=>$commit,expected=>$xv,expected_commit=>$xc,expected_at=>$xl,m32=>$m32,owner=>$e->[2],match=>("$v" eq $xv)?1:0} }
# every spelling of a MAX value that a row or tool must not carry: decimal, _-grouped, 2^k, 2^k-1, 2^a-2^b, (1 << a) - (1 << b)
sub pm_lits { my @out; for my $s (@_){ my $b=Math::BigInt->new($s); next if $b->is_zero; my $u="$b"; 1 while $u=~s/^(\d+)(\d{3})/$1_$2/; push @out,"$b"; push @out,$u if $u ne "$b";
    my $tz=0; my $x=$b->copy; while($x->copy->bmod(2)->is_zero){ $x->bdiv(2); $tz++ }
    my $y=$x->copy->binc; my $k=0; my $z=$y->copy; while($z->bcmp(1)>0 && $z->copy->bmod(2)->is_zero){ $z->bdiv(2); $k++ } my $pow2 = $z->bcmp(1)==0;
    if($pow2 && $tz==0){ push @out,"2^$k-1","2^$k - 1" } if($pow2 && $tz>0){ my $a=$k+$tz; push @out,"2^$a-2^$tz","2^$a - 2^$tz","(1 << $a) - (1 << $tz)" }
    if($x->bcmp(1)==0){ push @out,"2^$tz" } }
  my %s; grep { !$s{$_}++ } @out }
1;
