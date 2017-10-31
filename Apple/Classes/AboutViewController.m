//
//  AboutViewController.m
//  Piologie
//
//  Created by Sebastian Wedeniwski on 14.08.10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import <mach/mach.h>
#import <mach/mach_host.h>
#import "AboutViewController.h"
#import "IPadHelper.h"

@implementation AboutViewController

@synthesize logoView;
@synthesize numberOfDigitsTableView;
@synthesize aboutTextView;
@synthesize versionLabel;

-(natural_t) getFreeMemory {
  vm_size_t pagesize;
  vm_statistics_data_t vm_stat;
  mach_port_t host_port = mach_host_self();
  mach_msg_type_number_t host_size = sizeof(vm_statistics_data_t) / sizeof(integer_t);
  host_page_size(host_port, &pagesize);
  if (host_statistics(host_port, HOST_VM_INFO, (host_info_t)&vm_stat, &host_size) != KERN_SUCCESS) {
    NSLog(@"Failed to fetch vm statistics");
    return 0;
  }
  // Stats in bytes
  return vm_stat.free_count * pagesize;
}

-(void) createData {
  // view available memory
  //double m = ([self getFreeMemory]/1024.0)/1024.0;
  // Get the bundle path
  NSString *bPath = [[NSBundle mainBundle] bundlePath];
  NSString *settingsPath = [bPath stringByAppendingPathComponent:@"Settings.bundle"];
  NSString *plistFile = [settingsPath stringByAppendingPathComponent:@"Root.plist"];
  NSDictionary *settingsDictionary = [NSDictionary dictionaryWithContentsOfFile:plistFile];
  NSArray *settingsArray = [settingsDictionary objectForKey:@"PreferenceSpecifiers"];
  for (NSDictionary *item in settingsArray) {
    NSString *key = [item objectForKey:@"Key"];
    if ([key isEqualToString:@"max_constants_digits"]) {
      numberOfDigitsTableTitlesArray = [[NSArray alloc] initWithArray:[item objectForKey:@"Titles"]];
      numberOfDigitsTableValuesArray = [[NSArray alloc] initWithArray:[item objectForKey:@"Values"]];
      numberOfDigitsTableSectionTitle = [[item objectForKey:@"Title"] retain];
    } else if ([key isEqualToString:@"version"]) {
      versionLabel.text = [versionLabel.text stringByReplacingOccurrencesOfString:@"<version>" withString:[item objectForKey:@"DefaultValue"]];
    } else if ([key isEqualToString:@"copyright"]) {
      aboutTextView.text = [aboutTextView.text stringByReplacingOccurrencesOfString:@"<copyright>" withString:[item objectForKey:@"DefaultValue"]];
    }
  }
  NSUserDefaults *userDefaults = [NSUserDefaults standardUserDefaults];
  NSInteger m = [userDefaults integerForKey:@"max_constants_digits"];
  for (NSUInteger i = 0; i < [numberOfDigitsTableValuesArray count]; ++i) {
    NSInteger value = [[numberOfDigitsTableValuesArray objectAtIndex:i] integerValue];
    if (value == m) {
      numberOfDigitsTableSectionValue = [[numberOfDigitsTableTitlesArray objectAtIndex:i] retain];
      break;
    }
  }
}

#pragma mark -
#pragma mark Memory management

-(id)initWithNibName:(NSString *)nibNameOrNil bundle:(NSBundle *)nibBundleOrNil {
  self = [super initWithNibName:nibNameOrNil bundle:nibBundleOrNil];
  return self;
}

-(void)didReceiveMemoryWarning {
  [super didReceiveMemoryWarning];
}

-(void)viewDidUnload {
  [super viewDidUnload];
  numberOfDigitsTableView = nil;
  aboutTextView = nil;
  versionLabel = nil;
}

- (void)dealloc {
  [super dealloc];
  [numberOfDigitsTableTitlesArray release];
  [numberOfDigitsTableValuesArray release];
  [numberOfDigitsTableSectionTitle release];
  [numberOfDigitsTableSectionValue release];
  [logoView release];
  [numberOfDigitsTableView release];
  [aboutTextView release];
  [versionLabel release];
}

-(IBAction)loadBackView:(id)sender {
  [self dismissModalViewControllerAnimated:YES];
}

// Implement viewDidLoad to do additional setup after loading the view, typically from a nib.
- (void)viewDidLoad {
  [self createData];
  [super viewDidLoad];
  numberOfDigitsTableView.backgroundColor = [UIColor clearColor];
  if ([IPadHelper isIPad]) logoView.image = [UIImage imageNamed:@"Piologie@2x.png"];
}

/*-(BOOL)shouldAutorotateToInterfaceOrientation:(UIInterfaceOrientation)toInterfaceOrientation {
  return NO;
}*/

#pragma mark -
#pragma mark Table view data source

- (NSInteger)numberOfSectionsInTableView:(UITableView *)tableView {
  return 1;
}

- (NSInteger)tableView:(UITableView *)tableView numberOfRowsInSection:(NSInteger)section {
  return (selectNumberOfDigits == FALSE)? 1 : [numberOfDigitsTableTitlesArray count];
}

- (NSString *)tableView:(UITableView *)tableView titleForHeaderInSection:(NSInteger)section {
  return numberOfDigitsTableSectionTitle;
}

- (UITableViewCell *)tableView:(UITableView *)tableView cellForRowAtIndexPath:(NSIndexPath *)indexPath {
  static NSString *CellIdentifier = @"Cell";
  UITableViewCell *cell = [tableView dequeueReusableCellWithIdentifier:CellIdentifier];
  if (cell == nil) {
    cell = [[[UITableViewCell alloc] initWithStyle:UITableViewCellStyleDefault reuseIdentifier:CellIdentifier] autorelease];
  }
  if (selectNumberOfDigits == FALSE) {
    [[cell textLabel] setText:numberOfDigitsTableSectionValue];
    cell.accessoryType = UITableViewCellAccessoryDisclosureIndicator;
  } else {
    NSString *s = [numberOfDigitsTableTitlesArray objectAtIndex:indexPath.row];
    [[cell textLabel] setText:s];
    if ([s isEqualToString:numberOfDigitsTableSectionValue]) {
      cell.accessoryType = UITableViewCellAccessoryCheckmark;
      scrollPositionRow = indexPath.row;
    } else {
      cell.accessoryType = UITableViewCellAccessoryNone;
    }
  }
  return cell;
}

/*- (void)tableView:(UITableView *)tableView viewDidLoad {
  [super viewDidLoad];
  if (selectNumberOfDigits == TRUE) {
    [numberOfDigitsTableView scrollToRowAtIndexPath:[NSIndexPath indexPathForRow:scrollPositionRow inSection:0] atScrollPosition:UITableViewScrollPositionMiddle animated:NO];
  }
}*/

- (void)alertView:(UIAlertView *)alertView clickedButtonAtIndex:(NSInteger)buttonIndex {
  NSUserDefaults *userDefaults = [NSUserDefaults standardUserDefaults];
  [userDefaults setObject:[alertView buttonTitleAtIndex:buttonIndex] forKey:@"max_precision_warning"];
}

#pragma mark -
#pragma mark Table view delegate

- (void)tableView:(UITableView *)tableView didSelectRowAtIndexPath:(NSIndexPath *)indexPath {
  if (selectNumberOfDigits == FALSE) {
    selectNumberOfDigits = TRUE;
  } else {
    selectNumberOfDigits = FALSE;
    [numberOfDigitsTableSectionValue release];
    numberOfDigitsTableSectionValue = [[NSString alloc] initWithString:[numberOfDigitsTableTitlesArray objectAtIndex:indexPath.row]];
    NSUserDefaults *userDefaults = [NSUserDefaults standardUserDefaults];
    for (NSUInteger i = 0; i < [numberOfDigitsTableTitlesArray count]; ++i) {
      NSString *value = [numberOfDigitsTableTitlesArray objectAtIndex:i];
      if ([value isEqualToString:numberOfDigitsTableSectionValue]) {
        if ([[numberOfDigitsTableValuesArray objectAtIndex:i] integerValue] >= MAX_PRECISION_WARNING) {
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
        [userDefaults setObject:[numberOfDigitsTableValuesArray objectAtIndex:i] forKey:@"max_constants_digits"];
        break;
      }
    }
  }
  [numberOfDigitsTableView reloadData];
}

@end
