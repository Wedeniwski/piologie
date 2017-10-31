//
//  PiologieAppDelegate.h
//  Piologie
//
//  Created by Sebastian Wedeniwski on 08.08.10.
//  Copyright __MyCompanyName__ 2010. All rights reserved.
//

#import <UIKit/UIKit.h>

@class PiologieViewController;

@interface PiologieAppDelegate : NSObject <UIApplicationDelegate> {
    UIWindow *window;
    PiologieViewController *viewController;
}

@property (nonatomic, retain) IBOutlet UIWindow *window;
@property (nonatomic, retain) IBOutlet PiologieViewController *viewController;

@end

