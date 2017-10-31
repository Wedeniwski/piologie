//
//  CalculationViewController.m
//  Piologie
//
//  Created by Sebastian Wedeniwski on 08.08.10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import "CalculationViewController.h"
#import "constants.h"
#import "IPadHelper.h"

@implementation CalculationViewController

@synthesize numberOfDigitsSlider;
@synthesize numberOfDigitsLabel;
@synthesize calculationTimeLabel;
@synthesize positionOfViewDigitsSlider;
@synthesize calculationViewDigitsLabel;
@synthesize calculationResultTextView;
@synthesize waitView;


-(void)initData {
  numberOfDigits = 1000;
  digitsPerView = [IPadHelper isIPad]? 7700 : 1000;
  digitsPerLine = [IPadHelper isIPad]? 110 : 40;
  computationTime1000 = 0.0;

  NSUserDefaults *userDefaults = [NSUserDefaults standardUserDefaults];
  numberOfDigitsSlider.maximumValue = [userDefaults integerForKey:@"max_constants_digits"];
  numberOfDigitsSlider.minimumValue = 20;
  numberOfDigitsSlider.value = numberOfDigits;

  positionOfViewDigitsSlider.maximumValue = MAX(0, numberOfDigits-digitsPerView);
  positionOfViewDigitsSlider.minimumValue = 0;
  positionOfViewDigitsSlider.value = 0;
}

-(void) setConstantName:(NSString *)cName {
  constantName = [cName retain];
  constantNameLastCalc = @"";
}

-(double) estimatedCalculationTime {
  if ([constantName isEqualToString:CONSTANT_PI]) {
    // Pi:    1000 - 0.001028s
    // Pi:  100000 - 1.71s
    // Pi:  500000 - 13.87s
    // Pi: 1000000 - 31s
    // Pi: 5000000 - 457.6s
    double c = computationTime1000/0.001028;
    if (numberOfDigits < 100000) return 1.71*c*numberOfDigits/100000.0;
    else if (numberOfDigits < 500000) return 13.87*c*numberOfDigits/500000.0;
    else if (numberOfDigits < 1000000) return 31.0*c*numberOfDigits/1000000.0;
    return 457.6*c*numberOfDigits/5000000.0;
  } else if ([constantName isEqualToString:CONSTANT_SQRT2]) {
    // Sqrt2:    1000 - 0.000627s
    // Sqrt2:  100000 - 0.89s
    // Sqrt2:  500000 - 8.48s
    // Sqrt2: 1000000 - 16.88s
    // Sqrt2: 5000000 - 248.36s
    double c = computationTime1000/0.000627;
    if (numberOfDigits < 100000) return 0.89*c*numberOfDigits/100000.0;
    else if (numberOfDigits < 500000) return 8.48*c*numberOfDigits/500000.0;
    else if (numberOfDigits < 1000000) return 16.88*c*numberOfDigits/1000000.0;
    return 248.36*c*numberOfDigits/5000000.0;
  } else if ([constantName isEqualToString:CONSTANT_EXP1]) {
    // Exp1:    1000 - 0.000911s
    // Exp1:  100000 - 0.81s
    // Exp1:  500000 - 7.21s
    // Exp1: 1000000 - 13.46s
    // Exp1: 5000000 - 222.64s
    double c = computationTime1000/0.000911;
    if (numberOfDigits < 100000) return 0.81*c*numberOfDigits/100000.0;
    else if (numberOfDigits < 500000) return 7.21*c*numberOfDigits/500000.0;
    else if (numberOfDigits < 1000000) return 13.46*c*numberOfDigits/1000000.0;
    return 222.64*c*numberOfDigits/5000000.0;
  } else if ([constantName isEqualToString:CONSTANT_ZETA3]) {
    // Zeta3:    1000 - 0.001761s
    // Zeta3:  100000 - 5.4s
    // Zeta3:  500000 - 42.54s
    // Zeta3: 1000000 - 101.62s
    double c = computationTime1000/0.001761;
    if (numberOfDigits < 100000) return 5.4*c*numberOfDigits/100000.0;
    else if (numberOfDigits < 500000) return 42.54*c*numberOfDigits/500000.0;
    return 101.62*c*numberOfDigits/1000000.0;
  } else if ([constantName isEqualToString:CONSTANT_GAMMA]) {
    // Gamma:    1000 - 0.026186s
    // Gamma:   50000 - 59.82s
    // Gamma:  100000 - 149.57s
    // Gamma:  500000 - 646.61s
    // Gamma: 1000000 - 1252.13s
    double c = computationTime1000/0.026186;
    if (numberOfDigits < 50000) return 59.82*c*numberOfDigits/100000.0;
    else if (numberOfDigits < 100000) return 149.57*c*numberOfDigits/100000.0;
    else if (numberOfDigits < 500000) return 646.61*c*numberOfDigits/500000.0;
    return 1252.13*c*numberOfDigits/1000000.0;
  } else if ([constantName isEqualToString:CONSTANT_LN2]) {
    // Log2:    1000 - 0.002706s
    // Log2:  100000 - 11.28s
    // Log2:  500000 - 82.81s
    // Log2: 1000000 - 184.41s
    double c = computationTime1000/0.002706;
    if (numberOfDigits < 100000) return 11.28*c*numberOfDigits/100000.0;
    else if (numberOfDigits < 500000) return 82.81*c*numberOfDigits/500000.0;
    return 184.41*c*numberOfDigits/1000000.0;
  }
  return 0.0;
}

-(IBAction)updateView:(id)sender {
  numberOfDigits = numberOfDigitsSlider.value;
  if (![constantName isEqualToString:constantNameLastCalc] || numberOfDigits != numberOfDigitsLastCalc) {
    positionOfViewDigitsSlider.value = 0;
    positionOfViewDigitsSlider.maximumValue = (numberOfDigits > digitsPerView)? (numberOfDigits+digitsPerLine-digitsPerView) : 0;
    NSString *s = [[NSString alloc] initWithFormat:@"%d digits of constant %s:", numberOfDigits, [constantName UTF8String]];
    numberOfDigitsLabel.text = s;
    [s release];
    double c = [self estimatedCalculationTime];
    if (c < 0.1) {
      calculationTimeLabel.text = @"Calculation time...";
    } else {
      NSString *s = [[NSString alloc] initWithFormat:@"Calculation time...about %.0fs", c];
      calculationTimeLabel.text = s;
      [s release];
    }
    calculationResultTextView.text = @"";
    calculationViewDigitsLabel.text = @"";
    positionOfViewDigitsSlider.hidden = YES;
    waitView.hidden = NO;
  }
}

-(NSString *)outputCalcNumber:(int)pos {
  char d[50];
  const char* c = outputCalcUnformattedNumber();
  NSMutableString *strFormatted = [[[NSMutableString alloc] initWithCapacity:2*digitsPerView] autorelease];
  [strFormatted appendFormat:@"        %c", *c];
  while (*++c != '.') [strFormatted appendFormat:@"%c", *c];
  [strFormatted appendString:@".\n"];
  int length = strlen(++c);
  if (pos > 0 && length > pos) {
    [strFormatted setString:@""];
    c += pos;
    length -= pos;
  }
  bool sgn = (pos >= 100000);
  int pos2 = pos%100000;
  for (int i = 1; length > 0 && i <= digitsPerView; i += digitsPerLine) {
    int k = pos2+i;
    if (k >= 100000) {
      k -= 100000;
      sgn = true;
    }
    if (sgn) [strFormatted appendFormat:@"\'%05d:", k];
    else [strFormatted appendFormat:@"%6d:", k];
    for (int j = 0; length > 0 && j < digitsPerLine;) {
      if (length >= 40 && j+40 <= digitsPerLine) {
        d[0] = ' '; d[1] = c[0]; d[2] = c[1]; d[3] = c[2]; d[4] = c[3]; d[5] = c[4]; d[6] = c[5]; d[7] = c[6]; d[8] = c[7]; d[9] = c[8]; d[10] = c[9];
        d[11] = ' '; d[12] = c[10]; d[13] = c[11]; d[14] = c[12]; d[15] = c[13]; d[16] = c[14]; d[17] = c[15]; d[18] = c[16]; d[19] = c[17]; d[20] = c[18]; d[21] = c[19];
        d[22] = ' '; d[23] = c[20]; d[24] = c[21]; d[25] = c[22]; d[26] = c[23]; d[27] = c[24]; d[28] = c[25]; d[29] = c[26]; d[30] = c[27]; d[31] = c[28]; d[32] = c[29];
        d[33] = ' '; d[34] = c[30]; d[35] = c[31]; d[36] = c[32]; d[37] = c[33]; d[38] = c[34]; d[39] = c[35]; d[40] = c[36]; d[41] = c[37]; d[42] = c[38]; d[43] = c[39];
        d[44] = '\0';
        j += 40; c += 40; length -= 40;
      } else if (length >= 30 && j+30 <= digitsPerLine) {
        d[0] = ' '; d[1] = c[0]; d[2] = c[1]; d[3] = c[2]; d[4] = c[3]; d[5] = c[4]; d[6] = c[5]; d[7] = c[6]; d[8] = c[7]; d[9] = c[8]; d[10] = c[9];
        d[11] = ' '; d[12] = c[10]; d[13] = c[11]; d[14] = c[12]; d[15] = c[13]; d[16] = c[14]; d[17] = c[15]; d[18] = c[16]; d[19] = c[17]; d[20] = c[18]; d[21] = c[19];
        d[22] = ' '; d[23] = c[20]; d[24] = c[21]; d[25] = c[22]; d[26] = c[23]; d[27] = c[24]; d[28] = c[25]; d[29] = c[26]; d[30] = c[27]; d[31] = c[28]; d[32] = c[29];
        d[33] = '\0';
        j += 30; c += 30; length -= 30;
      } else if (length >= 20 && j+20 <= digitsPerLine) {
        d[0] = ' '; d[1] = c[0]; d[2] = c[1]; d[3] = c[2]; d[4] = c[3]; d[5] = c[4]; d[6] = c[5]; d[7] = c[6]; d[8] = c[7]; d[9] = c[8]; d[10] = c[9];
        d[11] = ' '; d[12] = c[10]; d[13] = c[11]; d[14] = c[12]; d[15] = c[13]; d[16] = c[14]; d[17] = c[15]; d[18] = c[16]; d[19] = c[17]; d[20] = c[18]; d[21] = c[19];
        d[22] = '\0';
        j += 19; c += 20; length -= 20;
      } else if (length >= 10 && j+10 <= digitsPerLine) {
        d[0] = ' '; d[1] = c[0]; d[2] = c[1]; d[3] = c[2]; d[4] = c[3]; d[5] = c[4]; d[6] = c[5]; d[7] = c[6]; d[8] = c[7]; d[9] = c[8]; d[10] = c[9];
        d[11] = '\0';
        j += 10; c += 10; length -= 10;
      } else {
        d[0] = ' ';
        int n = 1;
        for (; length > 0 && j < digitsPerLine; ++j, ++n, --length) d[n] = c[j];
        d[n] = '\0';
        c += n;
      }
      [strFormatted appendFormat:@"%s", d];
    }
    [strFormatted appendString:@"\n"];
  }
  return strFormatted;
}

-(IBAction)updateViewDigits:(id)sender {
  unsigned int n = positionOfViewDigitsSlider.value;
  n -= n%digitsPerLine;
  positionOfViewDigitsSlider.value = n;
  int m = n+digitsPerView;
  if (m > numberOfDigits) m = numberOfDigits;
  NSString *s = [[NSString alloc] initWithFormat:@"View digits from %d to %d:", (n+1), m];
  calculationViewDigitsLabel.text = s;
  [s release];
  /*const char* c = outputCalcNumber(n, digitsPerLine, digitsPerView);
  s = [[NSString alloc] initWithFormat:@"%s", c];
  //s = [[NSString alloc] initWithUTF8String:c];
  calculationResultTextView.text = s;
  [s release];*/
  calculationResultTextView.text = [self outputCalcNumber:n];
}

-(IBAction) calculate:(id)sender {
  [self updateView:sender];
  numberOfDigitsSlider.enabled = NO;
  
  positionOfViewDigitsSlider.hidden = YES;
  waitView.hidden = NO;
  if (![constantName isEqualToString:constantNameLastCalc] || numberOfDigits != numberOfDigitsLastCalc) {
    if ([constantName isEqualToString:CONSTANT_PI]) {
      calcPi(numberOfDigits);
    } else if ([constantName isEqualToString:CONSTANT_SQRT2]) {
      calcSqrt2(numberOfDigits);
    } else if ([constantName isEqualToString:CONSTANT_EXP1]) {
      calcExp1(numberOfDigits);
    } else if ([constantName isEqualToString:CONSTANT_ZETA3]) {
      calcZeta3(numberOfDigits);
    } else if ([constantName isEqualToString:CONSTANT_GAMMA]) {
      calcGamma(numberOfDigits);
    } else if ([constantName isEqualToString:CONSTANT_LN2]) {
      calcLn2(numberOfDigits);
    }
    double d1 = outputCalcDuration();
    double d2 = outputConversionDuration();
    if (computationTime1000 == 0.0 /*numberOfDigits == 1000*/) {
      computationTime1000 = d1+d2;
    }
    NSString *s = [[NSString alloc] initWithFormat:@"Calculation time: %.2fs (+ conversion time: %.2fs)", d1, d2];
    calculationTimeLabel.text = s;
    [s release];
    
    [self updateViewDigits:sender];

    [constantNameLastCalc release];
    constantNameLastCalc = [constantName retain];
    numberOfDigitsLastCalc = numberOfDigits;
  }
  waitView.hidden = YES;
  if (positionOfViewDigitsSlider.maximumValue > 0) {
    positionOfViewDigitsSlider.hidden = NO;
  }  
  numberOfDigitsSlider.enabled = YES;
}

-(IBAction)loadBackView:(id)sender {
  [self dismissModalViewControllerAnimated:YES];
}

-(NSString*) getArchimedesConstantInformation {
  return @"<html><head><script>document.ontouchmove = function(event) { if (document.body.scrollHeight == document.body.clientHeight) event.preventDefault(); }</script></head><body>The Archimedes\' constant π (written pi) is defined by the area enclosed by a circle of radius 1 while its circumference is equal to 2π. More general, the constant π is the ratio of a circle's area to the square of its radius <p><center><img src=\"Circle_Area.png\"></center></p> as also the ratio of a circle's circumference to its diameter <p><center><img src=\"Circle_Circumference.png\"></center></p>which was proved by Archimedes of Syracuse (287-212 BC).</body></html>";
}

-(NSString*) getPythagorasConstantInformation {
  return @"<html><head><script>document.ontouchmove = function(event) { if (document.body.scrollHeight == document.body.clientHeight) event.preventDefault(); }</script></head><body>The phrase \"Pythagoras\' constant\" is used for the square root of 2 which is the positive real number that, when multiplied by itself, gives the number 2.<br/>Geometrically, the square root of 2 is the length of the diagonal of a unit square. This follows from the Pythagorean theorem,<br/>i.e. 1<sup>2</sup> + 1<sup>2</sup> = 2.<p><center><img src=\"Square_root_of_2_triangle.png\"></center></p>According to the Greek philosopher Aristotle (384-322 BC), it was the Pythagoreans around 430 BC who first demonstrated the irrationality of the diagonal of the unit square. It was probably the first number discovered which did not fit in all their systems based on integers and fractions of integers.</body></html>";
}

-(NSString*) getEulersNumberInformation {
  return @"<html><head><script>document.ontouchmove = function(event) { if (document.body.scrollHeight == document.body.clientHeight) event.preventDefault(); }</script></head><body>Euler's number <i>e</i> is defined by the value of the following expression:<p><center><img src=\"EulersNumber.png\"></center></p> It can also be written as a Taylor series:<p><center><img src=\"EulersNumber2.png\"></center></p></body></html>";
}

-(NSString*) getAperyConstantInformation {
  return @"<html><head><script>document.ontouchmove = function(event) { if (document.body.scrollHeight == document.body.clientHeight) event.preventDefault(); }</script></head><body>The Apéry's constant is defined by the formula<p><center><img src=\"Apery's constant.png\"></center></p> where ζ is the Riemann zeta function. This constant was named for Roger Apéry, who in 1978 proved it to be irrational.<p>The constant ζ(3) has a probabilistic interpretation: Given three random integers, the probability that no factor exceeding 1 divides them all is 1/ζ(3).</p></body></html>";
}

-(NSString*) getEulerMascheroniConstantInformation {
  return @"<html><head><script>document.ontouchmove = function(event) { if (document.body.scrollHeight == document.body.clientHeight) event.preventDefault(); }</script></head><body>The Euler-Mascheroni constant is denoted by the lowercase Greek letter gamma (γ). It is defined by the limit<p><center><table><tr><td>γ =</td><td><img src=\"EulerGamma.png\"></td></tr></table></center></p>This constant measures the amount by which the partial sums of the harmonic series differ from the natural logarithm.<p>This constant appears frequently in number theory, for example, if a large integer <i>n</i> is divided by each integer 1 <u>&lt;</u> <i>k</i> <u>&lt;</u> <i>n</i>, then the average fraction by which the quotient <i>n</i>/<i>k</i> falls short of the next integer is not 1/2, but γ.</p></body></html>";
}

-(NSString*) getNaturalLogarithmOf2Information {
  return @"<html><head><script>document.ontouchmove = function(event) { if (document.body.scrollHeight == document.body.clientHeight) event.preventDefault(); }</script></head><body>The natural logarithm (short ln) is the logarithm to the base <i>e</i>, where <i>e</i> is Euler's Number. The natural logarithm of 2 is a transcendental quantity that arises often in decay problems, especially when half-lives are being converted to decay constants.<p>This value is the area under the hyperbola f(<i>x</i>)=1/<i>x</i> with 1 <u>&lt;</u> <i>x</i> <u>&lt;</u> 2.</br><center><img src=\"Natural_logarithm_integral.png\"></center></p></body></html>";
}

-(IBAction)loadInfoView:(id)sender {
  InformationViewController *controller = [[InformationViewController alloc] initWithNibName:@"InformationViewController" bundle:nil];
  NSString *s = @"<html><head><script>document.ontouchmove = function(event) { if (document.body.scrollHeight == document.body.clientHeight) event.preventDefault(); }</script></head><body>About: </body></html>";//[[NSString alloc] initWithFormat:@"<html><body>About: %s</body></html>", [constantName UTF8String]];
  if ([constantName isEqualToString:CONSTANT_PI]) {
    s = [self getArchimedesConstantInformation];
  } else if ([constantName isEqualToString:CONSTANT_SQRT2]) {
    s = [self getPythagorasConstantInformation];
  } else if ([constantName isEqualToString:CONSTANT_EXP1]) {
    s = [self getEulersNumberInformation];
  } else if ([constantName isEqualToString:CONSTANT_ZETA3]) {
    s = [self getAperyConstantInformation];
  } else if ([constantName isEqualToString:CONSTANT_GAMMA]) {
    s = [self getEulerMascheroniConstantInformation];
  } else if ([constantName isEqualToString:CONSTANT_LN2]) {
    s = [self getNaturalLogarithmOf2Information];
  }
  [controller setInformationContent:s];
  controller.modalTransitionStyle = UIModalTransitionStyleFlipHorizontal;
  [self presentModalViewController:controller animated:YES];
  [controller release];
  //[s release];
}

-(void)willRotateToInterfaceOrientation:(UIInterfaceOrientation)toInterfaceOrientation duration:(NSTimeInterval)duration {
  [super willRotateToInterfaceOrientation:toInterfaceOrientation duration:duration];
  if (toInterfaceOrientation == UIInterfaceOrientationPortrait || toInterfaceOrientation == UIInterfaceOrientationPortraitUpsideDown) {
    digitsPerLine = [IPadHelper isIPad]? 110 : 40;
    digitsPerView = [IPadHelper isIPad]? 7700 : 1000;
  } else {
    digitsPerLine = [IPadHelper isIPad]? 150 : 60;
    digitsPerView = [IPadHelper isIPad]? 7200 : 1000;
  }
  [self updateViewDigits:0];
}

#pragma mark -
#pragma mark Memory management

-(id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil {
  self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil];
  return self;
}

-(void)dealloc {
  [numberOfDigitsSlider release];
  [numberOfDigitsLabel release];
  [calculationTimeLabel release];
  [positionOfViewDigitsSlider release];
  [calculationViewDigitsLabel release];
  [calculationResultTextView release];
  [waitView release];
  [constantName release];
  [constantNameLastCalc release];
  [super dealloc];
}

-(void)didReceiveMemoryWarning {
  [super didReceiveMemoryWarning];
}

-(void)viewDidUnload {
  [super viewDidUnload];
  numberOfDigitsSlider = nil;
  numberOfDigitsLabel = nil;
  calculationTimeLabel = nil;
  positionOfViewDigitsSlider = nil;
  calculationViewDigitsLabel = nil;
  calculationResultTextView = nil;
  waitView = nil;
}

// Implement viewDidLoad to do additional setup after loading the view, typically from a nib.
- (void)viewDidLoad {
  [self initData];
  [super viewDidLoad];
  [self calculate:0];
}

-(BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)toInterfaceOrientation {
  return YES;
}

@end
