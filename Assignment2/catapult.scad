module cuboid(x,y,z,rx,ry,rz,lx,ly,lz){
    translate([lx,ly,lz]) 
    rotate([rx,ry,rz])
    cube([x,y,z], true);
}
module wheel(h, r){
    rotate([45,0,0])
    rotate([0,90,0])
    difference() {
        cylinder(h, r, r, true);
        translate([0, -r/2, 0]) cylinder(2*h, r/4, r/4, true);
        translate([0, r/2, 0]) cylinder(2*h, r/4, r/4, true);
        translate([-r/2, 0, 0]) cylinder(2*h, r/4, r/4, true);
        translate([r/2, 0, 0]) cylinder(2*h, r/4, r/4, true);
        
    }
    translate([-h/2,0,0]) 
    rotate([0,90,0])
    cylinder(h/2, r/4, r/4, true);
}

module all_wheel(w, l, h, r) {
    translate([-h/2,1,0])wheel(h,r);
    translate([w+h/2+1,1,0]) mirror([1,0,0]) wheel(h,r);
    translate([-h/2,l-1,0]) wheel(h,r);
    translate([w+h/2+1,l-1,0]) mirror([1,0,0])wheel(h,r);
}

module weight(l) {
    translate([-4.5, 11.5, -4.5]) cube([13, 0.5, 13]);
    translate([-4.5, 7.5, -4.5]) cube([13, 0.5, 13]);
    translate([-4, 4, -4]) cube([12, 12, 12]);
    translate([-4, 0, 0]) cube([12, 4, 4]);
    translate([0, 0, 4]) cube([4, 15, 4]);
    translate([0, 0, -4]) cube([4, l/2, 4]);
}
module bucket(l, b_r) {
    weight(l);
    translate([2,l,4]) 
    difference() {
        union(){
        translate([-2,-l,-4])cube([4, l, 4]);
        sphere(b_r);
            
        }
        sphere(b_r-1);
        translate([0,0,10]) cube([20,20,20], true);     
        rotate([30,-30,0])cylinder(20, 1, 1, true);    
        rotate([-30,-30,0])cylinder(20, 1, 1, true);    
        rotate([-30,30,0])cylinder(20, 1, 1, true);    
        rotate([30,30,0])cylinder(20, 1, 1, true);
        rotate([0,30,0])cylinder(20, 1, 1, true);
        rotate([0,-30,0])cylinder(20, 1, 1, true);
        rotate([-30,0,0])cylinder(20, 1, 1, true);
        rotate([30,0,0])cylinder(20, 1, 1, true);
        cylinder(20, 1, 1, true);
    }
}

module gear(handle, h, r) {
   rotate([0,90,0])
    union(){
   rotate([0,0,45]) union(){
   rotate([90,0,0]) cylinder(handle, 0.5, 0.5, true);
   rotate([0,90,0]) cylinder(handle, 0.5, 0.5, true);
   }
   rotate([0,90,0]) cylinder(handle, 0.5, 0.5, true);
   rotate([90,0,0]) cylinder(handle, 0.5, 0.5, true);
   cylinder(h, r, r, true);
   }
}

module base(l, w) {
    d = sqrt(w * w + l * l * 0.63); 
    cuboid(2,l,3,0,0,0,1,l/2,0);
    cuboid(2,l,3,0,0,0,w,l/2,0);
    cuboid(w,3,3,0,0,0,w/2,l*0.9,0);
    cuboid(w,3,3,0,0,0,w/2,2,0);

    cuboid(d,2,3,0,0,asin(l * 0.8/d),0.5+w/2,l/2,0);
    cuboid(d,2,3,0,0,-asin(l * 0.8/d),0.5+w/2,l/2,0);
}

module frame(w, l, s, r, h) {
    cuboid(2,s*2.3,3,60,0,0,1,cos(60)*s*1.15+1,sin(60)*s*1.15);
    cuboid(2,sqrt(2)*r*0.9,3,-45,0,0,1,s/2+cos(45)*r,r*0.9/2);
    cuboid(2,4,h,0,0,0,1,s*0.9,h/2);
    
    cuboid(2,s*2.3,3,60,0,0,w,cos(60)*s*1.15+1,sin(60)*s*1.15);
    cuboid(2,sqrt(2)*r*0.9,3,-45,0,0,w,s/2+cos(45)*r,r*0.9/2);
    cuboid(2,4,h,0,0,0,w,s*0.9,h/2);
}

module main(b_width=26, b_length=40, wheel_r=5, wheel_h=2, arm_angle=0, arm_size=60, gear_handle_length=16, gear_h=3, gear_r=5.5, gear_connector_r=1.75) {
union(){
all_wheel(b_width, b_length, wheel_h, wheel_r);
split_point = sqrt(3) / (sqrt(3) + 2.7) * b_length;
rest = b_length - split_point;

base(b_length, b_width);
height = sqrt(3)*split_point*0.85;
frame(b_width, b_length, split_point, rest, height);
    
translate([b_width/2+0.5, split_point*1.15, height*1.3]) 
    rotate([60,0,0])
    cube([b_width + 2, 2, 1], true);
    
translate([b_width/2, split_point*0.9, height*0.9]) 
    rotate([0,90,0])
    cylinder(b_width, 2, 2, true);

translate([b_width/2, split_point * 0.9, height*0.9]) 
    scale([1,1,1]) 
    rotate([arm_angle,0,0])
    translate([-2,(60-arm_size)*0.3-20,0]) 
    bucket(arm_size, 8);

translate([b_width/2,b_length*0.9*0.75, b_length*0.9*0.25]) 
    rotate([0,90,0])
    cylinder(b_width, gear_connector_r, gear_connector_r, true);
    
translate([3,b_length*0.9*0.75, b_length*0.9*0.25]) scale([1,0.75,0.75])gear(gear_handle_length, gear_h, gear_r);
translate([b_width-2,b_length*0.9*0.75, b_length*0.9*0.25]) scale([1,0.75,0.75])gear(gear_handle_length, gear_h, gear_r);
}
}
main(b_width=26, b_length=40, wheel_r=5, wheel_h=2, arm_angle=10, arm_size=60, gear_handle_length=16, gear_h=3, gear_r=5.5, gear_connector_r=1.75);