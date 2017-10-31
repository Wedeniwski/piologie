//
//  PiologieViewController.m
//  Piologie
//
//  Created by Sebastian Wedeniwski on 08.08.10.
//  Copyright __MyCompanyName__ 2010. All rights reserved.
//

#import "PiologieViewController.h"
#import "CalculationViewController.h"
#import "FactorizationViewController.h"
#import "InformationViewController.h"
#import "AboutViewController.h"

#define FACTORIZATION @"Factorization"
#define ABOUT         @"About"

@implementation PiologieViewController

@synthesize theTableView;

-(NSDictionary *)initialDefaults {
  NSArray *keys = [[[NSArray alloc] initWithObjects:@"max_constants_digits", @"max_precision_warning", nil] autorelease];
  NSArray *values = [[[NSArray alloc] initWithObjects:@"10000", @"Yes", nil] autorelease];
  return [[[NSDictionary alloc] initWithObjects: values forKeys: keys] autorelease];
}

-(void)alertView:(UIAlertView *)alertView clickedButtonAtIndex:(NSInteger)buttonIndex {
  NSUserDefaults *userDefaults = [NSUserDefaults standardUserDefaults];
  [userDefaults setObject:[alertView buttonTitleAtIndex:buttonIndex] forKey:@"max_precision_warning"];
}

-(void) createData {
  NSUserDefaults *userDefaults = [NSUserDefaults standardUserDefaults];
  [userDefaults registerDefaults:[self initialDefaults]];
  // ToDo: remove double code
  if ([userDefaults integerForKey:@"max_constants_digits"] >= MAX_PRECISION_WARNING) {
    if (![[userDefaults stringForKey:@"max_precision_warning"] isEqualToString:@"No"]) {
      UIAlertView *continueDialog = [[UIAlertView alloc]
                                     initWithTitle:@"High Precision Selection"
                                     message:@"This high precision might not work on all devices. Do not warn me again."
                                     delegate:self
                                     cancelButtonTitle:@"No"
                                     otherButtonTitles:@"Yes", nil];
      [continueDialog show];
      [continueDialog release];
    }
  }
  
  NSArray *numberTheory = [[NSArray alloc] initWithObjects:@"Number Theory", FACTORIZATION, nil];
  NSArray *constants = [[NSArray alloc] initWithObjects:@"Constants Computation", CONSTANT_PI, CONSTANT_SQRT2, CONSTANT_EXP1, CONSTANT_ZETA3, CONSTANT_GAMMA, CONSTANT_LN2, nil];
  NSArray *information = [[NSArray alloc] initWithObjects:@"Information", DOCUMENTAION, ABOUT, nil];
  theTableData = [[NSArray alloc] initWithObjects: numberTheory, constants, information, nil];
  [information release];
  [constants release];
  [numberTheory release];
}

/*
// The designated initializer. Override to perform setup that is required before the view is loaded.
- (id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil {
    if ((self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil])) {
        // Custom initialization
    }
    return self;
}
*/

/*
// Implement loadView to create a view hierarchy programmatically, without using a nib.
- (void)loadView {
}
*/


// Implement viewDidLoad to do additional setup after loading the view, typically from a nib.
- (void)viewDidLoad {
  [self createData];
  [super viewDidLoad];
  //theTableView.scrollEnabled = NO;
}

-(BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)toInterfaceOrientation {
  return YES;
}

- (void)factorizationViewControllerDidFinish:(FactorizationViewController *)controller {
  [self dismissModalViewControllerAnimated:YES];
  // refresh to remove previous selection
  [theTableView reloadData];
}

- (void)calculationViewControllerDidFinish:(CalculationViewController *)controller {
  [self dismissModalViewControllerAnimated:YES];
  // refresh to remove previous selection
  [theTableView reloadData];
}

- (void)informationViewControllerDidFinish:(InformationViewController *)controller {
  [self dismissModalViewControllerAnimated:YES];
  // refresh to remove previous selection
  [theTableView reloadData];
}

- (void)aboutViewControllerDidFinish:(AboutViewController *)controller {
  [self dismissModalViewControllerAnimated:YES];
  // refresh to remove previous selection
  [theTableView reloadData];
}

#pragma mark -
#pragma mark Table view data source

- (NSInteger)numberOfSectionsInTableView:(UITableView *)tableView {
  // Return the number of sections.
  return [theTableData count];
}


- (NSInteger)tableView:(UITableView *)tableView numberOfRowsInSection:(NSInteger)section {
  // Return the number of rows in the section.
  NSArray *data = [theTableData objectAtIndex:section];
  return [data count]-1;
}

- (NSString *)tableView:(UITableView *)tableView titleForHeaderInSection:(NSInteger)section {
  NSArray *data = [theTableData objectAtIndex:section];
  return [data objectAtIndex:0];
}

// Customize the appearance of table view cells.
- (UITableViewCell *)tableView:(UITableView *)tableView cellForRowAtIndexPath:(NSIndexPath *)indexPath {
  static NSString *CellIdentifier = @"Cell";
  UITableViewCell *cell = [tableView dequeueReusableCellWithIdentifier:CellIdentifier];
  if (cell == nil) {
    cell = [[[UITableViewCell alloc] initWithStyle:UITableViewCellStyleDefault reuseIdentifier:CellIdentifier] autorelease];
  }
  NSArray *data = [theTableData objectAtIndex:indexPath.section];
  // Configure the cell...
  [[cell textLabel] setText:[data objectAtIndex:(indexPath.row+1)]];
  cell.accessoryType = UITableViewCellAccessoryDisclosureIndicator;
  return cell;
}

#pragma mark -
#pragma mark Table view delegate

-(void)tableView:(UITableView *)tableView didSelectRowAtIndexPath:(NSIndexPath *)indexPath {
  NSArray *data = [theTableData objectAtIndex:indexPath.section];
  NSString *constantName = [data objectAtIndex:(indexPath.row+1)];
  if ([constantName isEqualToString:FACTORIZATION]) {
    FactorizationViewController *controller = [[FactorizationViewController alloc] initWithNibName:@"FactorizationViewController" bundle:nil];
    controller.modalTransitionStyle = UIModalTransitionStyleCoverVertical;
    [self presentModalViewController:controller animated:YES];
    [controller release];
  } else if ([constantName isEqualToString:DOCUMENTAION]) {
    /*NSString *filePath = [[NSBundle mainBundle] pathForResource:@"piologie" ofType:@"pdf"];
    NSURL *documentUrl = [NSURL fileURLWithPath:filePath];
    MFDocumentManager *aDocManager = [[MFDocumentManager alloc] initWithFileUrl:documentUrl];
    DocumentViewController_Kiosk *aDocViewController = [[DocumentViewController_Kiosk alloc] initWithDocumentManager:aDocManager];
    [aDocViewController setDocumentId:@"piologie"];   // We use the filename as an ID. You can use whaterver you like, like the id entry in a database or the hash of the document.
    //
    //	In this example we use a navigation controller to present the document view controller but you can present it
    //	as a modal viewcontroller or just show a single PDF right from the beginning
    [self presentModalViewController:aDocViewController animated:YES]; 
    
    //[[self navigationController]pushViewController:aDocViewController animated:YES];
    
    [aDocViewController release];
    [aDocManager release];*/

    InformationViewController *controller = [[InformationViewController alloc] initWithNibName:@"InformationViewController" bundle:nil];
    [controller setInformationContent:DOCUMENTAION];
    controller.modalTransitionStyle = UIModalTransitionStyleCoverVertical;
    [self presentModalViewController:controller animated:YES];
    [controller release];
  } else if ([constantName isEqualToString:ABOUT]) {
    AboutViewController *controller = [[AboutViewController alloc] initWithNibName:@"AboutViewController" bundle:nil];
    controller.modalTransitionStyle = UIModalTransitionStyleCoverVertical;
    [self presentModalViewController:controller animated:YES];
    [controller release];
  } else {
    CalculationViewController *controller = [[CalculationViewController alloc] initWithNibName:@"CalculationViewController" bundle:nil];
    [controller setConstantName:constantName];
    controller.modalTransitionStyle = UIModalTransitionStyleCoverVertical;
    [self presentModalViewController:controller animated:YES];
    [controller release];
  }
  [tableView deselectRowAtIndexPath:indexPath animated:NO];
}


- (void)didReceiveMemoryWarning {
	// Releases the view if it doesn't have a superview.
    [super didReceiveMemoryWarning];
	
	// Release any cached data, images, etc that aren't in use.
}

-(void)viewDidUnload {
  theTableView = nil;
}

-(void)dealloc {
  [super dealloc];
  [theTableData release];
  [theTableView release];
}

@end
