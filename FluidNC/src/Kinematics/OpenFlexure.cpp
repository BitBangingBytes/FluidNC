#include "OpenFlexure.h"

#include "../Machine/MachineConfig.h"
#include "../Limits.h"  // ambiguousLimit()
#include "../Machine/Homing.h"

#include "../Protocol.h"  // protocol_execute_realtime

#include <cmath>

/*
  ==================== How it Works ====================================
  On a delta machine, Grbl axis units are in radians
  The kinematics converts the cartesian moves in gcode into
  the radians to move the arms. The Grbl motion planner never sees
  the actual cartesian values.

  To make the moves straight and smooth on a delta, the cartesian moves
  are broken into small segments where the non linearity will not be noticed.
  This is similar to how Grbl draws arcs.

  For mpos reporting, the motor position in steps is proportional to arm angles 
  in radians, which is then converted to cartesian via the forward kinematics 
  transform. Arm angle 0 means horizontal.

  Positive angles are below horizontal.

  The machine's Z zero point in the kinematics is parallel to the arm axes.
  The offset of the Z distance from the arm axes to the end effector joints 
  at arm angle zero will be printed at startup on the serial port.

  Feedrate in gcode is in the cartesian units. This must be converted to the
  angles. This is done by calculating the segment move distance and the angle 
  move distance and applying that ration to the feedrate. 

  FYI: http://forums.trossenrobotics.com/tutorials/introduction-129/delta-robot-kinematics-3276/
  Better: http://hypertriangle.com/~alex/delta-robot-tutorial/


Default configuration

kinematics:
  OpenFlexure:

  TODO 
   - Constrain the geometry values to realistic values.
   - Implement $MI for dynamixel motors.

*/

namespace Kinematics {

    // trigonometric constants to speed up calculations
    const float sqrt3  = 1.732050807;
    const float dtr    = M_PI / (float)180.0;  // degrees to radians
    const float sin120 = sqrt3 / 2.0;
    const float cos120 = -0.5;
    const float tan60  = sqrt3;
    const float sin30  = 0.5;
    const float tan30  = 1.0 / sqrt3;

    double x_fac;
    double y_fac;
    double z_fac;
    double ttd[3][3];
    double tfd[3][3];

    static float last_angle[MAX_N_AXIS]     = { 0.0 };  // A place to save the previous motor angles for distance/feed rate calcs
    static float last_cartesian[MAX_N_AXIS] = { 0.0 };  // A place to save the previous motor angles for distance/feed rate calcs

    void OpenFlexure::group(Configuration::HandlerBase& handler) {
        handler.item("max_z_mm", _max_z, -10.0, 10.0); 
        handler.item("Backlash_X", _backlash_x);
        handler.item("Backlash_Y", _backlash_y);
        handler.item("Backlash_Z", _backlash_z);
        handler.item("X_Scale", _scale_x);
        handler.item("Y_Scale", _scale_y);
        handler.item("Flex_A", _flex_a);
        handler.item("Flex_B", _flex_b);
        handler.item("Flex_H", _flex_h);
        handler.item("Camera_Angle", _camera_angle, 0.0, 360.0);
    }

    void OpenFlexure::init() {
        // print a startup message to show the kinematics are enabled. Print the offset for reference
        log_info("Kinematic system:" << name() << " soft_limits:" << _softLimits);

        auto axes   = config->_axes;
        auto n_axis = config->_axes->_numberAxis;

        // warn about axissofy limits
        for (int axis = 0; axis < n_axis; axis++) {
            if (axes->_axis[axis]->_softLimits) {
                log_warn(" All soft_limits configured in axes should be false");
                break;
            }
        }
        /* Python code from OpenFlexure SangaDeltaStage
        x_fac: float = -1 * np.multiply( np.divide(2, np.sqrt(3)), np.divide(self.flex_b, self.flex_h) )
        y_fac: float = -1 * np.divide(self.flex_b, self.flex_h)
        z_fac: float = np.multiply(np.divide(1, 3), np.divide(self.flex_b, self.flex_a))
        */
        log_debug("Flex A: " << _flex_a << ", Flex B: " << _flex_b << ", Flex H: " << _flex_h);
        x_fac = ((-1.0 * ((2.0 / sqrt(3)) * ((double)_flex_b / (double)_flex_h))) * _scale_x);
        y_fac = ((-1.0 * ((double)_flex_b / (double)_flex_h)) * _scale_y);
        z_fac = ((1.0 / 3.0) * ((double)_flex_b / (double)_flex_a));
        log_debug("FAC: (" << x_fac << "," << y_fac << "," << z_fac << ")");

        tfd[0][0] = (-1.0 * x_fac);
        tfd[0][1] = x_fac;
        tfd[0][2] = 0.0;
        tfd[1][0] = (0.5 * y_fac);
        tfd[1][1] = (0.5 * y_fac);
        tfd[1][2] = (-1.0 * y_fac);
        tfd[2][0] = z_fac;
        tfd[2][1] = z_fac;
        tfd[2][2] = z_fac;   

        // computes the inverse of a matrix m
        double determinant = tfd[0][0] * (tfd[1][1] * tfd[2][2] - tfd[2][1] * tfd[1][2]) -
                             tfd[0][1] * (tfd[1][0] * tfd[2][2] - tfd[1][2] * tfd[2][0]) +
                             tfd[0][2] * (tfd[1][0] * tfd[2][1] - tfd[1][1] * tfd[2][0]);

        double invdet = 1 / determinant;

        ttd[0][0] = (tfd[1][1] * tfd[2][2] - tfd[2][1] * tfd[1][2]) * invdet;
        ttd[0][1] = (tfd[0][2] * tfd[2][1] - tfd[0][1] * tfd[2][2]) * invdet;
        ttd[0][2] = (tfd[0][1] * tfd[1][2] - tfd[0][2] * tfd[1][1]) * invdet;
        ttd[1][0] = (tfd[1][2] * tfd[2][0] - tfd[1][0] * tfd[2][2]) * invdet;
        ttd[1][1] = (tfd[0][0] * tfd[2][2] - tfd[0][2] * tfd[2][0]) * invdet;
        ttd[1][2] = (tfd[1][0] * tfd[0][2] - tfd[0][0] * tfd[1][2]) * invdet;
        ttd[2][0] = (tfd[1][0] * tfd[2][1] - tfd[2][0] * tfd[1][1]) * invdet;
        ttd[2][1] = (tfd[2][0] * tfd[0][1] - tfd[0][0] * tfd[2][1]) * invdet;
        ttd[2][2] = (tfd[0][0] * tfd[1][1] - tfd[1][0] * tfd[0][1]) * invdet;

        log_debug("tfd: (" << tfd[0][0] << "," << tfd[0][1] << "," << tfd[0][2] << ")");
        log_debug("tfd: (" << tfd[1][0] << "," << tfd[1][1] << "," << tfd[1][2] << ")");
        log_debug("tfd: (" << tfd[2][0] << "," << tfd[2][1] << "," << tfd[2][2] << ")");
        log_debug("ttd: (" << ttd[0][0] << "," << ttd[0][1] << "," << ttd[0][2] << ")");
        log_debug("ttd: (" << ttd[1][0] << "," << ttd[1][1] << "," << ttd[1][2] << ")");
        log_debug("ttd: (" << ttd[2][0] << "," << ttd[2][1] << "," << ttd[2][2] << ")");

        init_position();
    }

    void OpenFlexure::init_position() {
        float angles[MAX_N_AXIS]    = { 0.0, 0.0, 0.0 };
        float cartesian[MAX_N_AXIS] = { 0.0, 0.0, 0.0 };
        // Calculate the Z offset at the arm zero angles ...
        // Z offset is the z distance from the motor axes to the end effector axes at zero angle
        motors_to_cartesian(cartesian, angles, 3);  // Sets the cartesian values
        log_info("  Z Offset:" << cartesian[Z_AXIS]);
    }

    bool OpenFlexure::invalid_line(float* cartesian) {
        if (!_softLimits)
            return false;

        float motors[MAX_N_AXIS] = { 0.0, 0.0, 0.0 };

        if (!transform_cartesian_to_motors(motors, cartesian)) {
            limit_error();
            return true;
        }

        return false;
    }

    // TO DO. This is not supported yet. Other levels of protection will prevent "damage"
    bool OpenFlexure::invalid_arc(
        float* target, plan_line_data_t* pl_data, float* position, float center[3], float radius, size_t caxes[3], bool is_clockwise_arc) {
        return false;
    }

    // copied from Cartesian . Needs to be optimized for parallel delta.
    void OpenFlexure::constrain_jog(float* target, plan_line_data_t* pl_data, float* position) {
        // log_debug("Jog Test: from(" << position[X_AXIS] << ")"
        //                             << " to(" << target[X_AXIS] << ")");
        if (!_softLimits)
            return;

        float motors[MAX_N_AXIS] = { 0.0, 0.0, 0.0 };

        // Temp fix
        // If the target is reachable do nothing
        if (transform_cartesian_to_motors(motors, target)) {
            return;
        } else {
            log_warn("Kinematics soft limit jog rejection");
            copyAxes(target, position);
        }

        // TO DO better idea
        // loop back from the target in increments of  kinematic_segment_len_mm unitl the position is valid.
        // constrain to that target.
    }

    bool OpenFlexure::cartesian_to_motors(float* target, plan_line_data_t* pl_data, float* position) {
        float dx, dy, dz;  // distances in each cartesian axis
        float motor_angles[3] = { 0 };

        float seg_target[3];                    // The target of the current segment
        float feed_rate  = pl_data->feed_rate;  // save original feed rate
        bool  show_error = true;                // shows error once

        bool calc_ok = true;
        
        // auto n_axis = config->_axes->_numberAxis;
        // log_debug("Target Values"); 
        // for (int i=0; i < n_axis; i++) {
        //     log_debug(target[i] << "," );
        //     if ( (n_axis + 1) < n_axis)
        //         log_debug(",");
        // }
        // log_debug("Position Values"); 
        // for (int i=0; i < n_axis; i++) {
        //     log_debug(position[i] << "," );
        //     if ( (n_axis + 1) < n_axis)
        //         log_debug(",");
        // }
        
        // Check to see if an A, B or C axis value was provided and if so, act only on those.
        // We have to track these values and do a comparison test. 
        // I can probably use 'position' instead of motor_raw_offset
        bool raw_move = false;        
        // calc_ok = transform_cartesian_to_motors(motor_angles, position);  // Make sure motor angles reflect current location

        // if ( fabs(target[A_AXIS] - position[A_AXIS]) >= 0.001 ) {
        //     motor_angles[0] += (target[A_AXIS] - position[A_AXIS]);
        //     raw_move = true;
        // }
        // if ( fabs(target[B_AXIS] - position[B_AXIS]) >= 0.001 ) {
        //     motor_angles[1] += (target[B_AXIS] - position[B_AXIS]);
        //     raw_move = true;
        // }
        // if ( fabs(target[C_AXIS] - position[C_AXIS]) >= 0.001 ) {
        //     motor_angles[2] += (target[C_AXIS] - position[C_AXIS]);
        //     raw_move = true;
        // }
        // if (raw_move) {
        //     if(!mc_move_motors(motor_angles, pl_data)) {
        //         return false;
        //     }
        //     return true;
        // }

        // Test Hash (Skip all this old Delta code below)
        // target[X_AXIS] = position[X_AXIS] + ((target[X_AXIS] - position[X_AXIS]) * (1.0 / _scale_x));
        // target[Y_AXIS] = position[Y_AXIS] + ((target[Y_AXIS] - position[Y_AXIS]) * (1.0 / _scale_y));
        // calc_ok = transform_cartesian_to_motors(motor_angles, target);
        // if (!mc_move_motors(motor_angles, pl_data)) {
        //     return false;
        // }
        // return true;
        //Test Hash


        if (target[Z_AXIS] > _max_z) {
            log_debug("Kinematics error. Target:" << target[Z_AXIS] << " exceeds max_z:" << _max_z);
            return false;
        }

        //log_debug("Target (" << target[0] << "," << target[1] << "," << target[2]);

        calc_ok = transform_cartesian_to_motors(last_angle, position);
        if (!calc_ok) {
            log_warn("Kinematics error. Start position error (" << position[0] << "," << position[1] << "," << position[2] << ")");
            return false;
        }

        // Check the destination to see if it is in work area
        calc_ok = transform_cartesian_to_motors(motor_angles, target);
        if (!calc_ok) {
            log_warn("Kinematics error. Target unreachable (" << target[0] << "," << target[1] << "," << target[2] << ")");
            return false;
        }

        position[X_AXIS] += gc_state.coord_offset[X_AXIS];
        position[Y_AXIS] += gc_state.coord_offset[Y_AXIS];
        position[Z_AXIS] += gc_state.coord_offset[Z_AXIS];

        // calculate cartesian move distance for each axis
        dx         = target[X_AXIS] - position[X_AXIS];
        dy         = target[Y_AXIS] - position[Y_AXIS];
        dz         = target[Z_AXIS] - position[Z_AXIS];
        float dist = sqrt((dx * dx) + (dy * dy) + (dz * dz));

        // determine the number of segments we need	... round up so there is at least 1 (except when dist is 0)
        uint32_t segment_count = ceil(dist / _kinematic_segment_len_mm);

        float segment_dist = dist / ((float)segment_count);  // distance of each segment...will be used for feedrate conversion

        for (uint32_t segment = 1; segment <= segment_count; segment++) {
            if (sys.abort) {
                return true;
            }
            //log_debug("Segment:" << segment << " of " << segment_count);
            // determine this segment's target
            seg_target[X_AXIS] = position[X_AXIS] + (dx / float(segment_count) * segment);
            seg_target[Y_AXIS] = position[Y_AXIS] + (dy / float(segment_count) * segment);
            seg_target[Z_AXIS] = position[Z_AXIS] + (dz / float(segment_count) * segment);

            //log_debug("Segment target (" << seg_target[0] << "," << seg_target[1] << "," << seg_target[2] << ")");

            // calculate the delta motor angles
            bool calc_ok = transform_cartesian_to_motors(motor_angles, seg_target);

            if (!calc_ok) {
                if (show_error) {
                    log_error("Kinematic error motors (" << motor_angles[0] << "," << motor_angles[1] << "," << motor_angles[2] << ")");
                    show_error = false;
                }
                return false;
            }
            if (pl_data->motion.rapidMotion) {
                pl_data->feed_rate = feed_rate;
            } else {  //This might be slowing down the feed rate too much during jogs commands, revisit
                float delta_distance = three_axis_dist(motor_angles, last_angle);
                pl_data->feed_rate   = (feed_rate * delta_distance / segment_dist);  
            }

            // mc_line() returns false if a jog is cancelled.
            // In that case we stop sending segments to the planner.
            if (!mc_move_motors(motor_angles, pl_data)) {
                return false;
            }

            // save angles for next distance calc
            // This is after mc_line() so that we do not update
            // last_angle if the segment was discarded.
            memcpy(last_angle, motor_angles, sizeof(motor_angles));
        }
        return true;
    }

    void OpenFlexure::motors_to_cartesian(float* cartesian, float* motors, int n_axis) {
        // log_debug("motors_to_cartesian: Motors(" << motors[0] << "," << motors[1] << "," << motors[2] << ")");
        double result[3] = {
            ((tfd[0][0] * (double) motors[Z_AXIS]) + (tfd[0][1] * (double) motors[Y_AXIS]) + (tfd[0][2] * (double) motors[X_AXIS])),
            ((tfd[1][0] * (double) motors[Z_AXIS]) + (tfd[1][1] * (double) motors[Y_AXIS]) + (tfd[1][2] * (double) motors[X_AXIS])),
            ((tfd[2][0] * (double) motors[Z_AXIS]) + (tfd[2][1] * (double) motors[Y_AXIS]) + (tfd[2][2] * (double) motors[X_AXIS]))
        };
        cartesian[X_AXIS] = (float) result[0];
        cartesian[Y_AXIS] = (float) result[1];
        cartesian[Z_AXIS] = (float) result[2];
        // log_debug("motors_to_cartesian: Cartesian(" << cartesian[X_AXIS] << "," << cartesian[Y_AXIS] << "," << cartesian[Z_AXIS] << ")");
    }

    bool OpenFlexure::kinematics_homing(AxisMask& axisMask) {
        // only servos use custom homing. Steppers use limit switches
        // if (!_use_servos)
        //     false;
        return false; // We don't use servos

        auto axes   = config->_axes;
        auto n_axis = axes->_numberAxis;

        config->_axes->set_disable(false);

        // TODO deal with non kinematic axes above Z
        for (int axis = 0; axis < 3; axis++) {
            //set_motor_steps(axis, mpos_to_steps(axes->_axis[axis]->_homing->_mpos, axis));
            int32_t steps = mpos_to_steps(_homing_mpos, axis);
            set_motor_steps(axis, steps);
        }
        protocol_disable_steppers();
        return true;  // signal main code that this handled all homing
    }

    void OpenFlexure::releaseMotors(AxisMask axisMask, MotorMask motors) {}

    bool OpenFlexure::transform_cartesian_to_motors(float* motors, float* cartesian) {
        motors[0] = motors[1] = motors[2] = 0;
        bool calc_ok                      = false;

        if (cartesian[Z_AXIS] > _max_z) {
            log_debug("Kinematics transform error. Target:" << cartesian[Z_AXIS] << " exceeds max_z:" << _max_z);
            return false;
        }
        // Matrix multiply cartesian by ttd (Transform to Delta) to arrive at motor positions
        double result[3] = {
            ((ttd[0][0] * (double) cartesian[X_AXIS]) + (ttd[0][1] * (double) cartesian[Y_AXIS]) + (ttd[0][2] * (double) cartesian[Z_AXIS])),
            ((ttd[1][0] * (double) cartesian[X_AXIS]) + (ttd[1][1] * (double) cartesian[Y_AXIS]) + (ttd[1][2] * (double) cartesian[Z_AXIS])),
            ((ttd[2][0] * (double) cartesian[X_AXIS]) + (ttd[2][1] * (double) cartesian[Y_AXIS]) + (ttd[2][2] * (double) cartesian[Z_AXIS]))
        };
        // log_debug("C2M: ttd: (" << ttd[0][0] << "," << ttd[0][1] << "," << ttd[0][2] << ")");
        // log_debug("C2M: ttd: (" << ttd[1][0] << "," << ttd[1][1] << "," << ttd[1][2] << ")");
        // log_debug("C2M: ttd: (" << ttd[2][0] << "," << ttd[2][1] << "," << ttd[2][2] << ")");

        motors[0] = (float) result[2];
        motors[1] = (float) result[1];
        motors[2] = (float) result[0];
        calc_ok = true;
        // log_debug("cartesian_to_motors: Cartesian (" << cartesian[X_AXIS] << "," << cartesian[Y_AXIS] << "," << cartesian[Z_AXIS] << ")");
        // log_debug("cartesian_to_motors: Motors (" << motors[X_AXIS] << "," << motors[Y_AXIS] << "," << motors[Z_AXIS] << ")");
        return calc_ok;
    }

    // Determine the unit distance between (2) 3D points
    float OpenFlexure::three_axis_dist(float* point1, float* point2) {
        return sqrt(((point1[0] - point2[0]) * (point1[0] - point2[0])) + ((point1[1] - point2[1]) * (point1[1] - point2[1])) +
                    ((point1[2] - point2[2]) * (point1[2] - point2[2])));
    }

    // Configuration registration
    namespace {
        KinematicsFactory::InstanceBuilder<OpenFlexure> registration("OpenFlexure");
    }
}
