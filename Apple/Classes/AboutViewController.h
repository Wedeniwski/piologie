//
//  AboutViewController.h
//  Piologie
//
//  Created by Sebastian Wedeniwski on 14.08.10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import <UIKit/UIKit.h>

#define MAX_PRECISION_WARNING 5000000

@interface AboutViewController : UIViewController {
  BOOL selectNumberOfDigits;
  NSArray *numberOfDigitsTableTitlesArray;
  NSArray *numberOfDigitsTableValuesArray;
  NSString *numberOfDigitsTableSectionTitle;
  NSString *numberOfDigitsTableSectionValue;
  NSUInteger scrollPositionRow;

  IBOutlet UIImageView *logoView;
  IBOutlet UITableView *numberOfDigitsTableView;
  IBOutlet UITextView *aboutTextView;
  IBOutlet UILabel *versionLabel;
}

-(id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil;

-(IBAction)loadBackView:(id)sender;

@property (retain, nonatomic) UIImageView *logoView;
@property (retain, nonatomic) UITableView *numberOfDigitsTableView;
@property (retain, nonatomic) UITextView *aboutTextView;
@property (retain, nonatomic) UILabel *versionLabel;

@end
