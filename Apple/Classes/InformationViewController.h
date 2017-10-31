//
//  InformationViewController.h
//  Piologie
//
//  Created by Sebastian Wedeniwski on 11.08.10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import <UIKit/UIKit.h>
#import "PDFScrollView.h"

#define DOCUMENTAION  @"Documentation (German)"

@interface InformationViewController : UIViewController {
  int currentPage;
  float pageWidth;
  NSString *informationContent;
  PDFScrollView *pdfView;
  PDFScrollView *oldPdfView;
  UISwipeGestureRecognizer *swipeLeftRecognizer;
  UISwipeGestureRecognizer *swipeRightRecognizer;

  IBOutlet UIView *informationView;
  IBOutlet UILabel *pageLabel;
  IBOutlet UISlider *pageSlider;
  IBOutlet UIButton *previousPage;
  IBOutlet UIButton *nextPage;
}

-(id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil;

-(void)setInformationContent:(NSString*)content;
-(void)handleSwipeFrom:(UISwipeGestureRecognizer *)recognizer;

-(IBAction)loadBackView:(id)sender;
-(IBAction)pageChanged:(id)sender;
-(IBAction)pageChanging:(id)sender;
-(IBAction)previousPage:(id)sender;
-(IBAction)nextPage:(id)sender;

@property (retain, nonatomic) UIView *informationView;
@property (retain, nonatomic) UILabel *pageLabel;
@property (retain, nonatomic) UISlider *pageSlider;
@property (retain, nonatomic) UIButton *previousPage;
@property (retain, nonatomic) UIButton *nextPage;

@end
