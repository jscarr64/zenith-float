# Harper 2026-10-30 (Jeff ruling 6:00 PM CT, FIXLIST ZF-M27): one gold per row for i686 and x86_64. Shared by tools/lint.pl (BAD) and the pmrev7 counter.
# pm_target_specific(@f) returns the broken rule names for one 17-field record:
#  class      the class names a target (ends in -x86_64 or -i686)
#  gold       oracle_i686 carries a gold: token that oracle_x86_64 does not carry (a per-target gold file)
#  pin        required_result or refusal_variant names one target alone (after removing "i686 and x86_64",
#             "x86_64 and i686" and "on both targets"), i.e. a result pinned for one target only
# pm_m27_candidate(@f): a zf row with a RoundingMode parameter whose cells include a directed mode
#  (FIXLIST:145 "all 6 directed modes": every mode except ToEven) at a p that is not a multiple of 64.
sub pm_target_specific { my @f=@_; my @r;
  push @r,'class' if $f[5]=~/-(?:x86_64|i686)$/;
  my %g; $g{$1}=1 while $f[12]=~/(gold:\S+?)(?=[\s,;)]|$)/g; my $bad=0; while($f[13]=~/(gold:\S+?)(?=[\s,;)]|$)/g){ $bad=1 unless $g{$1} } push @r,'gold' if $bad;
  my $t="$f[10] $f[11]"; $t=~s/\b(?:i686 and x86_64|x86_64 and i686)\b//g; $t=~s/\bon both targets\b//g; push @r,'pin' if $t=~/\b(?:x86_64|i686)\b/;
  @r }
sub pm_m27_candidate { my @f=@_; return 0 unless $f[0] eq 'zf'; return 0 if $f[8]=~/^n\/a/ || $f[8] eq 'ToEven';
  return 1 if $f[5]=~/^rm-(None|Up|Down|ToZero|FromZero|ToOdd)$/; return 1 if $f[8]=~/all 7 modes/; 0 }
1;
