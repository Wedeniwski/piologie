//
//  FactorizationViewController.h
//  Piologie
//
//  Created by Sebastian Wedeniwski on 09.08.10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import <UIKit/UIKit.h>

@interface FactorizationViewController : UIViewController <UIAlertViewDelegate, UITextFieldDelegate> {
  NSArray *items;

  IBOutlet UIToolbar *navigationTitle;
  IBOutlet UITextField *numberForFactorizationTextField;
  IBOutlet UILabel *calculationTimeLabel;
  IBOutlet UIWebView *calculationResultWebView;
  IBOutlet UIActivityIndicatorView *waitView;
}

-(id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil;

-(IBAction)loadBackView:(id)sender;
-(IBAction)clearView:(id)sender;
-(IBAction)updateView:(id)sender;
-(IBAction)calculate:(id)sender;

@property (retain, nonatomic) UIToolbar *navigationTitle;
@property (retain, nonatomic) UITextField *numberForFactorizationTextField;
@property (retain, nonatomic) UILabel *calculationTimeLabel;
@property (retain, nonatomic) UIWebView *calculationResultWebView;
@property (retain, nonatomic) UIActivityIndicatorView *waitView;

@end
