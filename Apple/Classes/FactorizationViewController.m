//
//  FactorizationViewController.m
//  Piologie
//
//  Created by Sebastian Wedeniwski on 09.08.10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import "FactorizationViewController.h"
#import "factoring.h"

// challenging factorizations:
// 2^67-1 = 147573952589676412927
// 2^128+1 = 340282366920938463463374607431768211457


@implementation FactorizationViewController

@synthesize navigationTitle;
@synthesize numberForFactorizationTextField;
@synthesize calculationTimeLabel;
@synthesize calculationResultWebView;
@synthesize waitView;

#pragma mark -
#pragma mark Memory management

-(id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil {
  self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil];
  if (self != nil) {
    items = nil;
  }
  return self;
}

-(void)dealloc {
  [items release];
  [navigationTitle release];
  [numberForFactorizationTextField release];
  [calculationTimeLabel release];
  [calculationResultWebView release];
  [waitView release];
  [super dealloc];
}

-(void)didReceiveMemoryWarning {
  [super didReceiveMemoryWarning];
}

-(void)viewDidUnload {
  [super viewDidUnload];
  navigationTitle = nil;
  numberForFactorizationTextField = nil;
  calculationTimeLabel = nil;
  calculationResultWebView = nil;
  waitView = nil;
}

#pragma mark -
#pragma mark Actions

-(IBAction)loadBackView:(id)sender {
  [self dismissModalViewControllerAnimated:YES];
}

-(IBAction)clearView:(id)sender {
  [calculationResultWebView loadHTMLString:@"" baseURL:nil];
  calculationTimeLabel.text = @"";
  waitView.hidden = YES;
}

-(IBAction)updateView:(id)sender {
  waitView.hidden = NO;
}

- (void)alertView:(UIAlertView *)alertView clickedButtonAtIndex:(NSInteger)buttonIndex {
	NSString *buttonTitle = [alertView buttonTitleAtIndex:buttonIndex];
	if ([buttonTitle isEqualToString:@"Yes"]) {
    const char* c = longFactoring([numberForFactorizationTextField.text UTF8String]);
    NSString *html = [[NSString alloc] initWithFormat:@"%s", c];
    [calculationResultWebView loadHTMLString:html baseURL:nil];
    [html release];
    double d = outputFactoringDuration();
    NSString *s = [[NSString alloc] initWithFormat:@"Factorization in %.2fs:", d];
    calculationTimeLabel.text = s;
    [s release];
	}
  waitView.hidden = YES;
  navigationTitle.items = items;
}

-(IBAction)calculate:(id)sender {
  [numberForFactorizationTextField endEditing:NO];//resignFirstResponder];
}

// Implement viewDidLoad to do additional setup after loading the view, typically from a nib.
- (void)viewDidLoad {
  [self clearView:0];
  [super viewDidLoad];
  items = [navigationTitle.items retain];
}

-(BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)toInterfaceOrientation {
  return YES;
}

-(void)viewFactorization:(NSString *)factorization {
  if (!isFactoringComplete()) {
    UIAlertView *continueDialog = [[UIAlertView alloc]
                                   initWithTitle:@"Continue Factorization?"
                                   message:@"The factorization might need much longer. Do you want to continue?"
                                   delegate:self
                                   cancelButtonTitle:@"No"
                                   otherButtonTitles:@"Yes", nil];
    [continueDialog show];
    [continueDialog release];
  }
  [calculationResultWebView loadHTMLString:factorization baseURL:nil];
  [factorization release];
  double d = outputFactoringDuration();
  NSString *s = [[NSString alloc] initWithFormat:@"Factorization in %.2fs:", d];
  calculationTimeLabel.text = s;
  [s release];
  if (isFactoringComplete()) {
    waitView.hidden = YES;
    navigationTitle.items = items;
  }
}

-(void)calculateFactorization {
  NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
  const char* c = factoring([numberForFactorizationTextField.text UTF8String]);
  NSString *html = [[NSString alloc] initWithFormat:@"%s", c];
  [self performSelectorOnMainThread:@selector(viewFactorization:) withObject:html waitUntilDone:NO];
  [pool release];
}

-(void)textFieldDidEndEditing:(UITextField *)textField {
  [calculationResultWebView loadHTMLString:@"" baseURL:nil];
  calculationTimeLabel.text = @"";
  waitView.hidden = NO;
  navigationTitle.items = nil;
  [self performSelectorInBackground:@selector(calculateFactorization) withObject:nil];
}

-(BOOL)textFieldShouldReturn:(UITextField *)textField {
  [numberForFactorizationTextField endEditing:NO];
  return YES;
}

@end
