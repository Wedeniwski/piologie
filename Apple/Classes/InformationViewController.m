//
//  InformationViewController.m
//  Piologie
//
//  Created by Sebastian Wedeniwski on 11.08.10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import "InformationViewController.h"
//#import <QuartzCore/QuartzCore.h>

#define NUMBER_OF_PAGES 489

@implementation InformationViewController

@synthesize informationView;
@synthesize pageLabel;
@synthesize pageSlider;
@synthesize previousPage, nextPage;

#pragma mark -
#pragma mark Memory management

-(id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil {
  self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil];
  if (self != nil) {
    informationContent = nil;
    pdfView = nil;
    oldPdfView = nil;
  }
  return self;
}

-(void)dealloc {
  [informationContent release];
  [oldPdfView release];
  [pdfView release];
  [swipeLeftRecognizer release];
  [swipeRightRecognizer release];
  [informationView release];
  [pageLabel release];
  [pageSlider release];
  [previousPage release];
  [nextPage release];
  [super dealloc];
}

-(void)didReceiveMemoryWarning {
  [super didReceiveMemoryWarning];
}

-(void)viewDidUnload {
  [super viewDidUnload];
  swipeLeftRecognizer = nil;
  swipeRightRecognizer = nil;
  informationView = nil;
  pageLabel = nil;
  pageSlider = nil;
  previousPage = nil;
  nextPage = nil;
}

-(void)setInformationContent:(NSString*)content {
  [informationContent release];
  informationContent = [content retain];
}

/*- (void)accelerometer:(UIAccelerometer *)accelerometer didAccelerate:(UIAcceleration *)acceleration {
  // Get the rotation in radians.
  CGFloat rotation = (atan2(acceleration.x, acceleration.y) + M_PI);
  // Dampen it a bit.
  float div = 48 / (2 * M_PI);
  rotation = (int)(rotation * div) / div;
  // Animate the movement a bit.
  [UIView beginAnimations:nil context:NULL];
  [UIView setAnimationDuration:.2];
  //informationWebView.layer.transform = CATransform3DMakeRotation(rotation, 0, 0, 1);
  [UIView commitAnimations];
}*/

-(void)handleSwipeFrom:(UISwipeGestureRecognizer *)recognizer {
  if (recognizer.direction == UISwipeGestureRecognizerDirectionLeft) {
    [self nextPage:nextPage];
  } else  if (recognizer.direction == UISwipeGestureRecognizerDirectionRight) {
    [self previousPage:previousPage];
  }
}

-(void)viewPage:(int)animationDirection {
  if (pdfView != nil) {
    if (oldPdfView != nil) {
      [oldPdfView removeFromSuperview];
      [oldPdfView release];
    }
    oldPdfView = pdfView;
  }
  currentPage = (int)pageSlider.value;
  pageLabel.text = [NSString stringWithFormat:@"%d / %d", currentPage, NUMBER_OF_PAGES];
  pdfView = [[PDFScrollView alloc] initWithFrame:informationView.bounds atPage:currentPage];
  pdfView.autoresizingMask = UIViewAutoresizingFlexibleLeftMargin | UIViewAutoresizingFlexibleWidth | UIViewAutoresizingFlexibleRightMargin | UIViewAutoresizingFlexibleHeight;
  if (animationDirection == 0) {
    [informationView addSubview:pdfView];
  } else if (animationDirection < 0) {
    pdfView.frame = CGRectMake(pdfView.frame.origin.x+pageWidth, pdfView.frame.origin.y, pdfView.frame.size.width, pdfView.frame.size.height);
    [informationView addSubview:pdfView];
    [informationView sendSubviewToBack:pdfView];
    [UIView beginAnimations:nil context:NULL];
    [UIView setAnimationDuration:0.5];
    oldPdfView.frame = CGRectMake(oldPdfView.frame.origin.x-pageWidth, oldPdfView.frame.origin.y, oldPdfView.frame.size.width, oldPdfView.frame.size.height);
    //[UIView setAnimationTransition:UIViewAnimationTransitionCurlDown forView:oldPdfView cache:NO];
    [UIView commitAnimations];
    pdfView.frame = CGRectMake(pdfView.frame.origin.x-pageWidth, pdfView.frame.origin.y, pdfView.frame.size.width, pdfView.frame.size.height);
  } else {
    pdfView.frame = CGRectMake(pdfView.frame.origin.x-pageWidth, pdfView.frame.origin.y, pdfView.frame.size.width, pdfView.frame.size.height);
    [informationView addSubview:pdfView];
    [informationView sendSubviewToBack:pdfView];
    [UIView beginAnimations:nil context:NULL];
    [UIView setAnimationDuration:0.5];
    oldPdfView.frame = CGRectMake(oldPdfView.frame.origin.x+pageWidth, oldPdfView.frame.origin.y, oldPdfView.frame.size.width, oldPdfView.frame.size.height);
    //[UIView setAnimationTransition:UIViewAnimationTransitionCurlUp forView:oldPdfView cache:NO];
    [UIView commitAnimations];
    pdfView.frame = CGRectMake(pdfView.frame.origin.x+pageWidth, pdfView.frame.origin.y, pdfView.frame.size.width, pdfView.frame.size.height);
  }
}

// Implement viewDidLoad to do additional setup after loading the view, typically from a nib.
-(void)viewDidLoad {
  [super viewDidLoad];
  CGRect r = [[UIScreen mainScreen] bounds];
  pageWidth = r.size.width;
  if ([informationContent isEqualToString:DOCUMENTAION]) {
    pageSlider.value = 1.0;
    pageSlider.minimumValue = 1.0;
    pageSlider.maximumValue = NUMBER_OF_PAGES;
    [self viewPage:0];
    swipeLeftRecognizer = [[UISwipeGestureRecognizer alloc] initWithTarget:self action:@selector(handleSwipeFrom:)];
    swipeLeftRecognizer.direction = UISwipeGestureRecognizerDirectionLeft;
    [self.view addGestureRecognizer:swipeLeftRecognizer];
    swipeRightRecognizer = [[UISwipeGestureRecognizer alloc] initWithTarget:self action:@selector(handleSwipeFrom:)];
    swipeRightRecognizer.direction = UISwipeGestureRecognizerDirectionRight;
    [self.view addGestureRecognizer:swipeRightRecognizer];
  } else {
    pageSlider.hidden = YES;
    pageLabel.hidden = YES;
    previousPage.hidden = YES;
    nextPage.hidden = YES;
    NSURL *baseURL = [NSURL fileURLWithPath:[[NSBundle mainBundle] bundlePath]];
    UIWebView *webView = [[UIWebView alloc] initWithFrame:informationView.bounds];
    webView.autoresizingMask = UIViewAutoresizingFlexibleLeftMargin | UIViewAutoresizingFlexibleWidth | UIViewAutoresizingFlexibleRightMargin | UIViewAutoresizingFlexibleTopMargin | UIViewAutoresizingFlexibleHeight | UIViewAutoresizingFlexibleBottomMargin;
    [webView loadHTMLString:informationContent baseURL:baseURL];
    [informationView addSubview:webView];
  }
}

-(BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)toInterfaceOrientation {
  return YES;
}

-(void)willRotateToInterfaceOrientation:(UIInterfaceOrientation)toInterfaceOrientation duration:(NSTimeInterval)duration {
  CGRect r = [[UIScreen mainScreen] bounds];
  if (toInterfaceOrientation == UIInterfaceOrientationPortrait || toInterfaceOrientation == UIInterfaceOrientationPortraitUpsideDown) {
    pageWidth = r.size.width;
  } else {
    pageWidth = r.size.height;
  }
}

-(void)didRotateFromInterfaceOrientation:(UIInterfaceOrientation)fromInterfaceOrientation {
  if (pdfView != nil) [self viewPage:0];
}

-(IBAction)loadBackView:(id)sender {
  [self dismissModalViewControllerAnimated:YES];
}

-(IBAction)pageChanged:(id)sender {
  [self viewPage:currentPage-(int)pageSlider.value];
}

-(IBAction)pageChanging:(id)sender {
  int page = (int)pageSlider.value;
  pageLabel.text = [NSString stringWithFormat:@"%d / %d", page, NUMBER_OF_PAGES];
}

-(IBAction)previousPage:(id)sender {
  int page = (int)pageSlider.value;
  if (page > pageSlider.minimumValue) {
    pageSlider.value = page-1;
    [self pageChanged:sender];
  }
}

-(IBAction)nextPage:(id)sender {
  int page = (int)pageSlider.value;
  if (page < pageSlider.maximumValue) {
    pageSlider.value = page+1;
    [self pageChanged:sender];
  }
}

@end
