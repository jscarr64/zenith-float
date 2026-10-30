# Harper 2026-10-30 (Jeff ruling 6:00 PM CT, FIXLIST ZF-M27): one gold per row for i686 and x86_64. Shared by tools/lint.pl (BAD) and the pmrev7 counter.
# pm_target_specific(@f) returns the broken rule names for one 17-field record:
#  class      the class names a target (ends in -x86_64 or -i686)
#  gold       oracle_i686 carries a gold: token that oracle_x86_64 does not carry (a per-target gold file)
#  pin        required_result or refusal_variant names one target alone (after removing "i686 and x86_64",
#             "x86_64 and i686" and "on both targets"), i.e. a result pinned for one target only
# pm_m27_candidate(@f): a zf row with a RoundingMode parameter whose cells include a directed mode
#  (FIXLIST:145 "all 6 directed modes": ToEven, Up, Down, ToZero, FromZero, ToOdd; every mode except None, target-identical since Jeff 2026-10-29) at a p that is not a multiple of 64.
sub pm_target_specific { my @f=@_; my @r;
  push @r,'class' if $f[5]=~/-(?:x86_64|i686)$/;
  my %g; $g{$1}=1 while $f[12]=~/(gold:\S+?)(?=[\s,;)]|$)/g; my $bad=0; while($f[13]=~/(gold:\S+?)(?=[\s,;)]|$)/g){ $bad=1 unless $g{$1} } push @r,'gold' if $bad;
  my $t="$f[10] $f[11]"; $t=~s/\b(?:i686 and x86_64|x86_64 and i686)\b//g; $t=~s/\bon both targets\b//g; $t=~s/\b(?:on |for )?every 64-bit target\b//g; push @r,'pin' if $t=~/\b(?:x86_64|i686)\b/;
  @r }
sub pm_m27_candidate { my @f=@_; return 0 unless $f[0] eq 'zf'; return 0 if $f[8]=~/^n\/a/ || $f[5] eq 'rm-None';
  return 1 if $f[5]=~/^rm-(ToEven|Up|Down|ToZero|FromZero|ToOdd)$/; return 1 if $f[8]=~/all 7 modes/; return 1 if $f[5] eq 'prec-2^30-1'; 0 }
# pmrev8 P1 (pop-logs/PRINCIPLES.md): pm_target_dependent(@f) names every way a row's expected outcome could differ by target.
# Text inside a "Red until ... lands: ..." sentence or an "at 7955836 ..." clause describes today's bug, not the expected outcome, and is stripped first.
sub pm_outcome_text { my $t=join ' ',@_; $t=~s/[Rr]ed until .*?(?:\.\s|\.$|\)|$)//g; $t=~s/at 7955836[^.;)]*//g; $t=~s/\b(?:i686 and x86_64|x86_64 and i686)\b//g; $t=~s/\bon both targets\b//g; $t=~s/\b(?:on |for )?every 64-bit target\b//g; $t }
# pmrev10 P6 (pop-logs/PRINCIPLES.md P6, Jeff 2026-10-30 6:31 PM CT; FIXLIST ZF-M31): 64-bit minimum. pm_p6(@lines) returns one rule name per
#  32-bit target or 32-bit word reference in any field or the header (i686/i386/i586, "32-bit", a 32-bit target triple, target_pointer_width = "32", oracle_i686).
sub pm_p6 { my @r; for my $t (@_){ push @r,$1 while $t=~/((?<![A-Za-z0-9])i[3-6]86\b|\b32-bit\b|\bthumbv\w+|\barmv7\w*|target_pointer_width\s*(?:=|!=)\s*"?32"?)/gi } @r }
sub pm_target_dependent { my @f=@_; my @r; my $o=pm_outcome_text($f[10],$f[11]);
  push @r,'refusal' if pm_outcome_text($f[11])=~/\b(?:x86_64|i686|32-bit|64-bit)\b/;
  push @r,'limit' if $o=~/MAX_PRECISION\s*(?:=|is)\s*2\^(?!30-1\b)/ || $o=~/\b(?:32|64)-bit\b[^.;]*MAX_PRECISION|MAX_PRECISION[^.;]*\b(?:32|64)-bit\b/;
  my $tb="$f[9] $f[15]"; $tb=~s/\b(?:run )?on both targets\b//g; push @r,'budget' if $tb=~/\b(?:x86_64|i686)\b|-only shard|pm-max-(?:i686|x86_64)/;
  push @r,'dependent' if $o=~/\b(?:32|64)-bit targets?\b|word size|usize::BITS|target_pointer_width|target-dependent|per-target|depends on the target/i;
  @r }
1;
