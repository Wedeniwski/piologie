//
//  PiologieViewController.h
//  Piologie
//
//  Created by Sebastian Wedeniwski on 08.08.10.
//  Copyright __MyCompanyName__ 2010. All rights reserved.
//

#import <UIKit/UIKit.h>

@interface PiologieViewController : UIViewController {
  NSArray *theTableData;
  IBOutlet UITableView *theTableView;
}

@property (retain, nonatomic) UITableView *theTableView;

@end

