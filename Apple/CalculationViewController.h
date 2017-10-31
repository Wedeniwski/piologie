//
//  CalculationViewController.h
//  Piologie
//
//  Created by Sebastian Wedeniwski on 08.08.10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import <UIKit/UIKit.h>
#import "InformationViewController.h"

#define CONSTANT_PI     @"Archimedes\' Constant Pi"
#define CONSTANT_SQRT2  @"Pythagoras\' Constant"
#define CONSTANT_EXP1   @"Euler\'s Number e"
#define CONSTANT_ZETA3  @"Apery\'s Constant"
#define CONSTANT_GAMMA  @"Euler-Mascheroni Constant"
#define CONSTANT_LN2    @"Natural Logarithm of 2"

@interface CalculationViewController : UIViewController {
  IBOutlet UISlider *numberOfDigitsSlider;
  IBOutlet UILabel *numberOfDigitsLabel;
  IBOutlet UILabel *calculationTimeLabel;
  IBOutlet UISlider *positionOfViewDigitsSlider;
  IBOutlet UILabel *calculationViewDigitsLabel;
  IBOutlet UITextView *calculationResultTextView;
  IBOutlet UIActivityIndicatorView *waitView;
  unsigned int numberOfDigits;
  unsigned int numberOfDigitsLastCalc;
  unsigned int digitsPerLine;
  unsigned int digitsPerView;
  double computationTime1000;
  NSString *constantName;
  NSString *constantNameLastCalc;
}

-(id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil;

-(void)setConstantName:(NSString *)cName;

-(IBAction)updateView:(id)sender;
-(IBAction)updateViewDigits:(id)sender;
-(IBAction)calculate:(id)sender;
-(IBAction)loadBackView:(id)sender;
-(IBAction)loadInfoView:(id)sender;

@property (retain, nonatomic) UISlider *numberOfDigitsSlider;
@property (retain, nonatomic) UILabel *numberOfDigitsLabel;
@property (retain, nonatomic) UILabel *calculationTimeLabel;
@property (retain, nonatomic) UISlider *positionOfViewDigitsSlider;
@property (retain, nonatomic) UILabel *calculationViewDigitsLabel;
@property (retain, nonatomic) UITextView *calculationResultTextView;
@property (retain, nonatomic) UIActivityIndicatorView *waitView;

@end
