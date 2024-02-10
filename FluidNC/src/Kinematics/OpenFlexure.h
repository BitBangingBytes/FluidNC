// Copyright (c) 2021 -	Bart Dring
// Use of this source code is governed by a GPLv3 license that can be found in the LICENSE file.

#pragma once

/*

	This implements Parallel Delta Kinematics

*/

#include "Kinematics.h"
#include "Cartesian.h"

// M_PI is not defined in standard C/C++ but some compilers
// support it anyway.  The following suppresses Intellisense
// problem reports.
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

namespace Kinematics {

    class OpenFlexure : public Cartesian {
    public:
        OpenFlexure() = default;

        OpenFlexure(const OpenFlexure&)            = delete;
        OpenFlexure(OpenFlexure&&)                 = delete;
        OpenFlexure& operator=(const OpenFlexure&) = delete;
        OpenFlexure& operator=(OpenFlexure&&)      = delete;

        // Kinematic Interface
        virtual void init() override;
        virtual void init_position() override;
        //bool canHome(AxisMask& axisMask) override;
        bool cartesian_to_motors(float* target, plan_line_data_t* pl_data, float* position) override;
        void motors_to_cartesian(float* cartesian, float* motors, int n_axis) override;
        bool transform_cartesian_to_motors(float* motors, float* cartesian) override;
        //bool soft_limit_error_exists(float* cartesian) override;
        bool         kinematics_homing(AxisMask& axisMask) override;
        virtual void constrain_jog(float* cartesian, plan_line_data_t* pl_data, float* position) override;
        virtual bool invalid_line(float* cartesian) override;
        virtual bool invalid_arc(float*            target,
                                 plan_line_data_t* pl_data,
                                 float*            position,
                                 float             center[3],
                                 float             radius,
                                 size_t            caxes[3],
                                 bool              is_clockwise_arc) override;

        void releaseMotors(AxisMask axisMask, MotorMask motors) override;

        // Configuration handlers:
        //void         validate() const override {}
        virtual void group(Configuration::HandlerBase& handler) override;
        void         afterParse() override {}

        // Name of the configurable. Must match the name registered in the cpp file.
        virtual const char* name() const override { return "OpenFlexure"; }

        ~OpenFlexure() {}

    private:
        //  Config items Using geometry names from the published kinematics rather than typical Fluid Style
        // To make the math easier to compare with the code
        // float of_rf = 70.0;  // The length of the crank arm on the motor
        // float of_f  = 179.437;
        // float of_re = 133.50;
        // float of_e  = 86.603;

        float _kinematic_segment_len_mm = 1.0;  // the maximun segment length the move is broken into
        bool  _softLimits               = false;
        float _homing_mpos              = 0.0;
        float _max_z                    = 0.0;
        // bool  _use_servos               = true;  // servo use a special homing
        int _backlash_x                 = 0;
        int _backlash_y                 = 0;
        int _backlash_z                 = 0;
        float _scale_x                  = 1.0;
        float _scale_y                  = 1.0;
        int _flex_a                     = 50;
        int _flex_b                     = 50;
        int _flex_h                     = 80;
        float _camera_angle               = 0.0;

        // bool  delta_calcAngleYZ(float x0, float y0, float z0, float& theta);
        float three_axis_dist(float* point1, float* point2);

        double ttd[3][3];
        double tfd[3][3];
        double x_fac = 0.0;
        double y_fac = 0.0;
        double z_fac = 0.0;

    protected:
    };
}  //  namespace Kinematics
