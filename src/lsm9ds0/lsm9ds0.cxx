/*
 * Author: Jon Trulson <jtrulson@ics.com>
 * Copyright (c) 2015 Intel Corporation.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <unistd.h>
#include <iostream>
#include <stdexcept>
#include <string.h>

// Madgwick
#include <math.h>
// Definitions
#define sampleFreqDef   512.0f          // sample frequency in Hz
#define betaDef         0.1f            // 2 * proportional gain

#include "lsm9ds0.hpp"

using namespace upm;
using namespace std;

// AHRS algorithm update
Madgwick::Madgwick() {
  beta = betaDef;
  q0 = 1.0f;
  q1 = 0.0f;
  q2 = 0.0f;
  q3 = 0.0f;
  invSampleFreq = 1.0f / sampleFreqDef;
  anglesComputed = 0;
}

void Madgwick::update(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz) {
  float recipNorm;
  float s0, s1, s2, s3;
  float qDot1, qDot2, qDot3, qDot4;
  float hx, hy;
  float _2q0mx, _2q0my, _2q0mz, _2q1mx, _2bx, _2bz, _4bx, _4bz, _2q0, _2q1, _2q2, _2q3, _2q0q2, _2q2q3, q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;

  // Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
  if((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f)) {
    updateIMU(gx, gy, gz, ax, ay, az);
    return;
  }

  // Convert gyroscope degrees/sec to radians/sec
  gx *= 0.0174533f;
  gy *= 0.0174533f;
  gz *= 0.0174533f;

  // Rate of change of quaternion from gyroscope
  qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
  qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy);
  qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx);
  qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx);

  // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
  if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

    // Normalise accelerometer measurement
    recipNorm = invSqrt(ax * ax + ay * ay + az * az);
    ax *= recipNorm;
    ay *= recipNorm;
    az *= recipNorm;

    // Normalise magnetometer measurement
    recipNorm = invSqrt(mx * mx + my * my + mz * mz);
    mx *= recipNorm;
    my *= recipNorm;
    mz *= recipNorm;

    // Auxiliary variables to avoid repeated arithmetic
    _2q0mx = 2.0f * q0 * mx;
    _2q0my = 2.0f * q0 * my;
    _2q0mz = 2.0f * q0 * mz;
    _2q1mx = 2.0f * q1 * mx;
    _2q0 = 2.0f * q0;
    _2q1 = 2.0f * q1;
    _2q2 = 2.0f * q2;
    _2q3 = 2.0f * q3;
    _2q0q2 = 2.0f * q0 * q2;
    _2q2q3 = 2.0f * q2 * q3;
    q0q0 = q0 * q0;
    q0q1 = q0 * q1;
    q0q2 = q0 * q2;
    q0q3 = q0 * q3;
    q1q1 = q1 * q1;
    q1q2 = q1 * q2;
    q1q3 = q1 * q3;
    q2q2 = q2 * q2;
    q2q3 = q2 * q3;
    q3q3 = q3 * q3;

    // Reference direction of Earth's magnetic field
    hx = mx * q0q0 - _2q0my * q3 + _2q0mz * q2 + mx * q1q1 + _2q1 * my * q2 + _2q1 * mz * q3 - mx * q2q2 - mx * q3q3;
    hy = _2q0mx * q3 + my * q0q0 - _2q0mz * q1 + _2q1mx * q2 - my * q1q1 + my * q2q2 + _2q2 * mz * q3 - my * q3q3;
    _2bx = sqrtf(hx * hx + hy * hy);
    _2bz = -_2q0mx * q2 + _2q0my * q1 + mz * q0q0 + _2q1mx * q3 - mz * q1q1 + _2q2 * my * q3 - mz * q2q2 + mz * q3q3;
    _4bx = 2.0f * _2bx;
    _4bz = 2.0f * _2bz;

    // Gradient decent algorithm corrective step
    s0 = -_2q2 * (2.0f * q1q3 - _2q0q2 - ax) + _2q1 * (2.0f * q0q1 + _2q2q3 - ay) - _2bz * q2 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q3 + _2bz * q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
    s1 = _2q3 * (2.0f * q1q3 - _2q0q2 - ax) + _2q0 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q1 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + _2bz * q3 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q2 + _2bz * q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q3 - _4bz * q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
    s2 = -_2q0 * (2.0f * q1q3 - _2q0q2 - ax) + _2q3 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q2 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + (-_4bx * q2 - _2bz * q0) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q1 + _2bz * q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q0 - _4bz * q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
    s3 = _2q1 * (2.0f * q1q3 - _2q0q2 - ax) + _2q2 * (2.0f * q0q1 + _2q2q3 - ay) + (-_4bx * q3 + _2bz * q1) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q0 + _2bz * q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
    recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
    s0 *= recipNorm;
    s1 *= recipNorm;
    s2 *= recipNorm;
    s3 *= recipNorm;

    // Apply feedback step
    qDot1 -= beta * s0;
    qDot2 -= beta * s1;
    qDot3 -= beta * s2;
    qDot4 -= beta * s3;
  }

  // Integrate rate of change of quaternion to yield quaternion
  q0 += qDot1 * invSampleFreq;
  q1 += qDot2 * invSampleFreq;
  q2 += qDot3 * invSampleFreq;
  q3 += qDot4 * invSampleFreq;

  // Normalise quaternion
  recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
  q0 *= recipNorm;
  q1 *= recipNorm;
  q2 *= recipNorm;
  q3 *= recipNorm;
  anglesComputed = 0;
}

//-------------------------------------------------------------------------------------------
// IMU algorithm update

void Madgwick::updateIMU(float gx, float gy, float gz, float ax, float ay, float az) {
  float recipNorm;
  float s0, s1, s2, s3;
  float qDot1, qDot2, qDot3, qDot4;
  float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3;

  // Convert gyroscope degrees/sec to radians/sec
  gx *= 0.0174533f;
  gy *= 0.0174533f;
  gz *= 0.0174533f;

  // Rate of change of quaternion from gyroscope
  qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
  qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy);
  qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx);
  qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx);

  // Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
  if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

    // Normalise accelerometer measurement
    recipNorm = invSqrt(ax * ax + ay * ay + az * az);
    ax *= recipNorm;
    ay *= recipNorm;
    az *= recipNorm;

    // Auxiliary variables to avoid repeated arithmetic
    _2q0 = 2.0f * q0;
    _2q1 = 2.0f * q1;
    _2q2 = 2.0f * q2;
    _2q3 = 2.0f * q3;
    _4q0 = 4.0f * q0;
    _4q1 = 4.0f * q1;
    _4q2 = 4.0f * q2;
    _8q1 = 8.0f * q1;
    _8q2 = 8.0f * q2;
    q0q0 = q0 * q0;
    q1q1 = q1 * q1;
    q2q2 = q2 * q2;
    q3q3 = q3 * q3;

    // Gradient decent algorithm corrective step
    s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
    s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
    s2 = 4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
    s3 = 4.0f * q1q1 * q3 - _2q1 * ax + 4.0f * q2q2 * q3 - _2q2 * ay;
    recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
    s0 *= recipNorm;
    s1 *= recipNorm;
    s2 *= recipNorm;
    s3 *= recipNorm;

    // Apply feedback step
    qDot1 -= beta * s0;
    qDot2 -= beta * s1;
    qDot3 -= beta * s2;
    qDot4 -= beta * s3;
  }

  // Integrate rate of change of quaternion to yield quaternion
  q0 += qDot1 * invSampleFreq;
  q1 += qDot2 * invSampleFreq;
  q2 += qDot3 * invSampleFreq;
  q3 += qDot4 * invSampleFreq;

  // Normalise quaternion
  recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
  q0 *= recipNorm;
  q1 *= recipNorm;
  q2 *= recipNorm;
  q3 *= recipNorm;
  anglesComputed = 0;
}

//-------------------------------------------------------------------------------------------
// Fast inverse square-root
// See: http://en.wikipedia.org/wiki/Fast_inverse_square_root

float Madgwick::invSqrt(float x) {
  float halfx = 0.5f * x;
  float y = x;
  long i = *(long*)&y;
  i = 0x5f3759df - (i>>1);
  y = *(float*)&i;
  y = y * (1.5f - (halfx * y * y));
  y = y * (1.5f - (halfx * y * y));
  return y;
}

//-------------------------------------------------------------------------------------------

void Madgwick::computeAngles()
{
  roll = atan2f(q0*q1 + q2*q3, 0.5f - q1*q1 - q2*q2);
  pitch = asinf(-2.0f * (q1*q3 - q0*q2));
  yaw = atan2f(q1*q2 + q0*q3, 0.5f - q2*q2 - q3*q3);
  anglesComputed = 1;
}


LSM9DS0::LSM9DS0(int bus, bool raw, uint8_t gAddress, uint8_t xmAddress) :
  m_i2cG(bus, raw), m_i2cXM(bus, raw), m_gpioG_INT(0), m_gpioG_DRDY(0),
  m_gpioXM_GEN1(0), m_gpioXM_GEN2(0)
{
  m_gAddr = gAddress;
  m_xmAddr = xmAddress;

  m_accelX = 0.0;
  m_accelY = 0.0;
  m_accelZ = 0.0;
  
  m_gyroX = 0.0;
  m_gyroY = 0.0;
  m_gyroZ = 0.0;
  
  m_magX = 0.0;
  m_magY = 0.0;
  m_magZ = 0.0;
  
  m_temp = 0.0;

  m_accelScale = 0.0;
  m_gyroScale = 0.0;
  m_magScale = 0.0;

  beta = betaDef;
  q0 = 1.0f;
  q1 = 0.0f;
  q2 = 0.0f;
  q3 = 0.0f;
  invSampleFreq = 1.0f / sampleFreqDef;
  anglesComputed = 0;

  mraa::Result rv;
  if ( (rv = m_i2cG.address(m_gAddr)) != mraa::SUCCESS)
    {
      throw std::runtime_error(string(__FUNCTION__) +
                               ": Could not initialize Gyro i2c address");
      return;
    }

  if ( (rv = m_i2cXM.address(m_xmAddr)) != mraa::SUCCESS)
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Could not initialize XM i2c address");
      return;
    }
}

LSM9DS0::~LSM9DS0()
{
  uninstallISR(INTERRUPT_G_INT);
  uninstallISR(INTERRUPT_G_DRDY);
  uninstallISR(INTERRUPT_XM_GEN1);
  uninstallISR(INTERRUPT_XM_GEN2);
}

bool LSM9DS0::init()
{
  // Init the gyroscope

  // power up
  if (!setGyroscopePowerDown(false))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to wake up gyro");
      return false;
    }
  
  // enable all axes
  if (!setGyroscopeEnableAxes(CTRL_REG1_G_YEN |CTRL_REG1_G_XEN |
                              CTRL_REG1_G_ZEN))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to enable gyro axes");
      return false;
    }
  
  // set gyro ODR
  if (!setGyroscopeODR(G_ODR_95_25))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to set gyro ODR");
      return false;
    }

  // set gyro scale
  if (!setGyroscopeScale(G_FS_245))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to set gyro scale");
      return false;
    }

  // Init the accelerometer

  // power up and set ODR
  if (!setAccelerometerODR(XM_AODR_100))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to set accel ODR");
      return false;
    }

  // enable all axes
  if (!setAccelerometerEnableAxes(CTRL_REG1_XM_AXEN |CTRL_REG1_XM_AYEN |
                                  CTRL_REG1_XM_AZEN))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to enable accel axes");
      return false;
    }
  
  // set scaling rate
  if (!setAccelerometerScale(XM_AFS_2))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to set accel scale");
      return false;
    }
  
  // temperature sensor

  // enable the temperature sensor
  if (!enableTemperatureSensor(true))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to enable temp sensor");
      return false;
    }

  // Init the magnetometer
  
  // set mode (this also powers it up if not XM_MD_POWERDOWN)
  if (!setMagnetometerMode(XM_MD_CONTINUOUS))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to set mag scale");
      return false;
    }

  // turn LPM off
  if (!setMagnetometerLPM(false))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to disable mag LPM");
      return false;
    }

  // set resolution
  if (!setMagnetometerResolution(XM_RES_LOW))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to set mag res");
      return false;
    }
  
  // set ODR
  if (!setMagnetometerODR(XM_ODR_12_5))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to set mag ODR");
      return false;
    }
  
  // set scale
  if (!setMagnetometerScale(XM_MFS_2))
    {
      throw std::runtime_error(string(__FUNCTION__) + 
                               ": Unable to set mag scale");
      return false;
    }

  return true;
}


void LSM9DS0::update()
{
  updateGyroscope();
  updateAccelerometer();
  updateMagnetometer();
  updateTemperature();
}

void LSM9DS0::updateGyroscope()
{
  uint8_t buffer[6];

  memset(buffer, 0, 6);
  readRegs(DEV_GYRO, REG_OUT_X_L_G, buffer, 6);

  int16_t x, y, z;

  x =  ( (buffer[1] << 8) | buffer[0] );
  y =  ( (buffer[3] << 8) | buffer[2] );
  z =  ( (buffer[5] << 8) | buffer[4] );

  m_gyroX = float(x);
  m_gyroY = float(y);
  m_gyroZ = float(z);
}

void LSM9DS0::updateAccelerometer()
{
  uint8_t buffer[6];

  memset(buffer, 0, 6);
  readRegs(DEV_XM, REG_OUT_X_L_A, buffer, 6);

  int16_t x, y, z;

  x =  ( (buffer[1] << 8) | buffer[0] );
  y =  ( (buffer[3] << 8) | buffer[2] );
  z =  ( (buffer[5] << 8) | buffer[4] );

  m_accelX = float(x);
  m_accelY = float(y);
  m_accelZ = float(z);
}

void LSM9DS0::updateMagnetometer()
{
  uint8_t buffer[6];

  memset(buffer, 0, 6);
  readRegs(DEV_XM, REG_OUT_X_L_M, buffer, 6);

  int16_t x, y, z;

  x =  ( (buffer[1] << 8) | buffer[0] );
  y =  ( (buffer[3] << 8) | buffer[2] );
  z =  ( (buffer[5] << 8) | buffer[4] );

  m_magX = float(x);
  m_magY = float(y);
  m_magZ = float(z);
}

void LSM9DS0::updateTemperature()
{
  uint8_t buffer[2];

  memset(buffer, 0, 2);
  readRegs(DEV_XM, REG_OUT_TEMP_L_XM, buffer, 2);

  //  cerr << "HIGH: " << int(buffer[1]) << " LOW: " << int(buffer[0]) << endl;

  // 12b signed
  int16_t temp = ( (buffer[1] << 8) | (buffer[0] ) );
  if (temp & 0x0800)
    {
      temp &= ~0x0800;
      temp *= -1;
    }

  m_temp = float(temp);
}

uint8_t LSM9DS0::readReg(DEVICE_T dev, uint8_t reg)
{
  mraa::I2c *device;

  switch(dev)
    {
    case DEV_GYRO: device = &m_i2cG; break;
    case DEV_XM:   device = &m_i2cXM; break;
    default:
      throw std::logic_error(string(__FUNCTION__) + 
                             ": Internal error, invalid device specified");
      return 0;
    }

  return device->readReg(reg);
}

void LSM9DS0::readRegs(DEVICE_T dev, uint8_t reg, uint8_t *buffer, int len)
{
  mraa::I2c *device;

  switch(dev)
    {
    case DEV_GYRO: device = &m_i2cG; break;
    case DEV_XM:   device = &m_i2cXM; break;
    default:
      throw std::logic_error(string(__FUNCTION__) + 
                             ": Internal error, invalid device specified");
      return;
    }

  // We need to set the high bit of the register to enable
  // auto-increment mode for reading multiple registers in one go.
  device->readBytesReg(reg | m_autoIncrementMode, buffer, len);
}

bool LSM9DS0::writeReg(DEVICE_T dev, uint8_t reg, uint8_t val)
{
  mraa::I2c *device;

  switch(dev)
    {
    case DEV_GYRO: device = &m_i2cG; break;
    case DEV_XM:   device = &m_i2cXM; break;
    default:
      throw std::logic_error(string(__FUNCTION__) + 
                             ": Internal error, invalid device specified");
      return false;
    }

  mraa::Result rv;
  if ((rv = device->writeReg(reg, val)) != mraa::SUCCESS)
    {
      throw std::runtime_error(std::string(__FUNCTION__) +
                               ": I2c.writeReg() failed");
      return false;
    } 
  
  return true;
}

bool LSM9DS0::setGyroscopePowerDown(bool enable)
{
  uint8_t reg = readReg(DEV_GYRO, REG_CTRL_REG1_G);

  if (enable)
    reg &= ~CTRL_REG1_G_PD;
  else
    reg |= CTRL_REG1_G_PD;

  return writeReg(DEV_GYRO, REG_CTRL_REG1_G, reg);
}

bool LSM9DS0::setGyroscopeEnableAxes(uint8_t axes)
{
  uint8_t reg = readReg(DEV_GYRO, REG_CTRL_REG1_G);

  // filter out any non-axis related data from arg
  axes &= (CTRL_REG1_G_YEN | CTRL_REG1_G_XEN |  CTRL_REG1_G_ZEN);

  // clear them in the register
  reg &= ~(CTRL_REG1_G_YEN | CTRL_REG1_G_XEN |  CTRL_REG1_G_ZEN);

  // now add them
  reg |= axes;

  return writeReg(DEV_GYRO, REG_CTRL_REG1_G, reg);
}

bool LSM9DS0::setGyroscopeODR(G_ODR_T odr)
{
  uint8_t reg = readReg(DEV_GYRO, REG_CTRL_REG1_G);

  reg &= ~(_CTRL_REG1_G_ODR_MASK << _CTRL_REG1_G_ODR_SHIFT);

  reg |= (odr << _CTRL_REG1_G_ODR_SHIFT);
  
  return writeReg(DEV_GYRO, REG_CTRL_REG1_G, reg);
}

bool LSM9DS0::setGyroscopeScale(G_FS_T scale)
{
  uint8_t reg = readReg(DEV_GYRO, REG_CTRL_REG4_G);

  reg &= ~(_CTRL_REG4_G_FS_MASK << _CTRL_REG4_G_FS_SHIFT);

  reg |= (scale << _CTRL_REG4_G_FS_SHIFT);

  if (!writeReg(DEV_GYRO, REG_CTRL_REG4_G, reg))
    {
      return false;
    }

  // store scaling factor (mDeg/s/LSB)

  switch (scale)
    {
    case G_FS_245:
      m_gyroScale = 8.75;
      break;

    case G_FS_500:
      m_gyroScale = 17.5;
      break;

    case G_FS_2000:
      m_gyroScale = 70.0;
      break;

    default: // should never occur, but...
      m_gyroScale = 0.0;        // set a safe, though incorrect value
      throw std::logic_error(string(__FUNCTION__) + 
                             ": internal error, unsupported scale");
      break;
    }

  return true;
}

bool LSM9DS0::setAccelerometerEnableAxes(uint8_t axes)
{
  uint8_t reg = readReg(DEV_XM, REG_CTRL_REG1_XM);

  // filter out any non-axis related data from arg
  axes &= (CTRL_REG1_XM_AXEN | CTRL_REG1_XM_AYEN | CTRL_REG1_XM_AZEN);

  // clear them in the register
  reg &= ~(CTRL_REG1_XM_AXEN | CTRL_REG1_XM_AYEN | CTRL_REG1_XM_AZEN);

  // now add them
  reg |= axes;

  return writeReg(DEV_XM, REG_CTRL_REG1_XM, reg);
}

bool LSM9DS0::setAccelerometerODR(XM_AODR_T odr)
{
  uint8_t reg = readReg(DEV_XM, REG_CTRL_REG1_XM);

  reg &= ~(_CTRL_REG1_XM_AODR_MASK << _CTRL_REG1_XM_AODR_SHIFT);

  reg |= (odr << _CTRL_REG1_XM_AODR_SHIFT);
  
  return writeReg(DEV_XM, REG_CTRL_REG1_XM, reg);
}

bool LSM9DS0::setAccelerometerScale(XM_AFS_T scale)
{
  uint8_t reg = readReg(DEV_XM, REG_CTRL_REG2_XM);

  reg &= ~(_CTRL_REG2_XM_AFS_MASK << _CTRL_REG2_XM_AFS_SHIFT);

  reg |= (scale << _CTRL_REG2_XM_AFS_SHIFT);

  if (!writeReg(DEV_XM, REG_CTRL_REG2_XM, reg))
    {
      return false;
    }

  // store scaling factor
  
  switch (scale)
    {
    case XM_AFS_2:
      m_accelScale = 0.061; 
      break;

    case XM_AFS_4:
      m_accelScale = 0.122 ;
      break;

    case XM_AFS_6:
      m_accelScale = 0.183 ;
      break;

    case XM_AFS_8:
      m_accelScale = 0.244 ;
      break;

    case XM_AFS_16:
      m_accelScale = 0.732 ;
      break;

    default: // should never occur, but...
      m_accelScale = 0.0;        // set a safe, though incorrect value
      throw std::logic_error(string(__FUNCTION__) + 
                             ": internal error, unsupported scale");
      break;
    }

  return true;
}

bool LSM9DS0::setMagnetometerResolution(XM_RES_T res)
{
  uint8_t reg = readReg(DEV_XM, REG_CTRL_REG5_XM);

  reg &= ~(_CTRL_REG5_XM_RES_MASK << _CTRL_REG5_XM_RES_SHIFT);

  reg |= (res << _CTRL_REG5_XM_RES_SHIFT);
  
  return writeReg(DEV_XM, REG_CTRL_REG5_XM, reg);
}

bool LSM9DS0::setMagnetometerODR(XM_ODR_T odr)
{
  uint8_t reg = readReg(DEV_XM, REG_CTRL_REG5_XM);

  reg &= ~(_CTRL_REG5_XM_ODR_MASK << _CTRL_REG5_XM_ODR_SHIFT);

  reg |= (odr << _CTRL_REG5_XM_ODR_SHIFT);
  
  return writeReg(DEV_XM, REG_CTRL_REG5_XM, reg);
}

bool LSM9DS0::setMagnetometerMode(XM_MD_T mode)
{
  uint8_t reg = readReg(DEV_XM, REG_CTRL_REG7_XM);

  reg &= ~(_CTRL_REG7_XM_MD_MASK << _CTRL_REG7_XM_MD_SHIFT);

  reg |= (mode << _CTRL_REG7_XM_MD_SHIFT);
  
  return writeReg(DEV_XM, REG_CTRL_REG7_XM, reg);
}

bool LSM9DS0::setMagnetometerLPM(bool enable)
{
  uint8_t reg = readReg(DEV_XM, REG_CTRL_REG7_XM);

  if (enable)
    reg |= CTRL_REG7_XM_MLP;
  else
    reg &= ~CTRL_REG7_XM_MLP;
  
  return writeReg(DEV_XM, REG_CTRL_REG7_XM, reg);
}

bool LSM9DS0::setMagnetometerScale(XM_MFS_T scale)
{
  uint8_t reg = readReg(DEV_XM, REG_CTRL_REG6_XM);

  reg &= ~(_CTRL_REG6_XM_MFS_MASK << _CTRL_REG6_XM_MFS_SHIFT);

  reg |= (scale << _CTRL_REG6_XM_MFS_SHIFT);

  if (!writeReg(DEV_XM, REG_CTRL_REG6_XM, reg))
    {
      return false;
    }

  // store scaling factor
  
  switch (scale)
    {
    case XM_MFS_2:
      m_magScale = 0.08;
      break;

    case XM_MFS_4:
      m_magScale = 0.16;
      break;

    case XM_MFS_8:
      m_magScale = 0.32;
      break;

    case XM_MFS_12:
      m_magScale = 0.48;
      break;

    default: // should never occur, but...
      m_magScale = 0.0;        // set a safe, though incorrect value
      throw std::logic_error(string(__FUNCTION__) + 
                             ": internal error, unsupported scale");
      break;
    }

  return true;
}

void LSM9DS0::getAccelerometer(float *x, float *y, float *z)
{
  if (x)
    *x = (m_accelX * m_accelScale) / 1000.0;

  if (y)
    *y = (m_accelY * m_accelScale) / 1000.0;

  if (z)
    *z = (m_accelZ * m_accelScale) / 1000.0;
}

void LSM9DS0::getGyroscope(float *x, float *y, float *z)
{
  if (x)
    *x = (m_gyroX * m_gyroScale) / 1000.0;

  if (y)
    *y = (m_gyroY * m_gyroScale) / 1000.0;

  if (z)
    *z = (m_gyroZ * m_gyroScale) / 1000.0;
}

void LSM9DS0::getMagnetometer(float *x, float *y, float *z)
{
  if (x)
    *x = (m_magX * m_magScale) / 1000.0;

  if (y)
    *y = (m_magY * m_magScale) / 1000.0;

  if (z)
    *z = (m_magZ * m_magScale) / 1000.0;
}

#ifdef JAVACALLBACK
float *LSM9DS0::getAccelerometer()
{
  float *v = new float[3];
  getAccelerometer(&v[0], &v[1], &v[2]);
  return v;
}

float *LSM9DS0::getGyroscope()
{
  float *v = new float[3];
  getGyroscope(&v[0], &v[1], &v[2]);
  return v;
}

float *LSM9DS0::getMagnetometer()
{
  float *v = new float[3];
  getMagnetometer(&v[0], &v[1], &v[2]);
  return v;
}
#endif

float LSM9DS0::getTemperature()
{
  // This might be wrong... The datasheet does not provide enough info
  // to calculate the temperature given a specific sensor reading.  So
  // - with 12b resolution, signed, and 8 degrees/per LSB, we come up
  // with the following.  Then scale up and we get a number that seems
  // pretty close.
  return (((m_temp / 2048.0) * 8.0) * 100.0);
}

bool LSM9DS0::enableTemperatureSensor(bool enable)
{
  uint8_t reg = readReg(DEV_XM, REG_CTRL_REG5_XM);

  if (enable)
    reg |= CTRL_REG5_XM_TEMP_EN;
  else
    reg &= ~CTRL_REG5_XM_TEMP_EN;

  return writeReg(DEV_XM, REG_CTRL_REG5_XM, reg);
}

uint8_t LSM9DS0::getGyroscopeStatus()
{
  return readReg(DEV_GYRO, REG_STATUS_REG_G);
}

uint8_t LSM9DS0::getMagnetometerStatus()
{
  return readReg(DEV_XM, REG_STATUS_REG_M);
}

uint8_t LSM9DS0::getAccelerometerStatus()
{
  return readReg(DEV_XM, REG_STATUS_REG_A);
}

uint8_t LSM9DS0::getGyroscopeInterruptConfig()
{
  return readReg(DEV_GYRO, REG_INT1_CFG_G);
}

bool LSM9DS0::setGyroscopeInterruptConfig(uint8_t enables)
{
  return writeReg(DEV_GYRO, REG_INT1_CFG_G, enables);
}

uint8_t LSM9DS0::getGyroscopeInterruptSrc()
{
  return readReg(DEV_GYRO, REG_INT1_SRC_G);
}

uint8_t LSM9DS0::getMagnetometerInterruptControl()
{
  return readReg(DEV_XM, REG_INT_CTRL_REG_M);
}

bool LSM9DS0::setMagnetometerInterruptControl(uint8_t enables)
{
  return writeReg(DEV_XM, REG_INT_CTRL_REG_M, enables);
}

uint8_t LSM9DS0::getMagnetometerInterruptSrc()
{
  return readReg(DEV_XM, REG_INT_SRC_REG_M);
}

uint8_t LSM9DS0::getInterruptGen1()
{
  return readReg(DEV_XM, REG_INT_GEN_1_REG);
}

bool LSM9DS0::setInterruptGen1(uint8_t enables)
{
  return writeReg(DEV_XM, REG_INT_GEN_1_REG, enables);
}

uint8_t LSM9DS0::getInterruptGen1Src()
{
  return readReg(DEV_XM, REG_INT_GEN_1_SRC);
}

uint8_t LSM9DS0::getInterruptGen2()
{
  return readReg(DEV_XM, REG_INT_GEN_2_REG);
}

bool LSM9DS0::setInterruptGen2(uint8_t enables)
{
  return writeReg(DEV_XM, REG_INT_GEN_2_REG, enables);
}

uint8_t LSM9DS0::getInterruptGen2Src()
{
  return readReg(DEV_XM, REG_INT_GEN_2_SRC);
}

#if defined(SWIGJAVA) || defined (JAVACALLBACK)
void LSM9DS0::installISR(INTERRUPT_PINS_T intr, int gpio, mraa::Edge level,
			 jobject runnable)
{
  // delete any existing ISR and GPIO context
  uninstallISR(intr);

  // greate gpio context
  getPin(intr) = new mraa::Gpio(gpio);

  getPin(intr)->dir(mraa::DIR_IN);
  getPin(intr)->isr(level, runnable);

}
#else
void LSM9DS0::installISR(INTERRUPT_PINS_T intr, int gpio, mraa::Edge level, 
                         void (*isr)(void *), void *arg)
{
  // delete any existing ISR and GPIO context
  uninstallISR(intr);

  // greate gpio context
  getPin(intr) = new mraa::Gpio(gpio);

  getPin(intr)->dir(mraa::DIR_IN);
  getPin(intr)->isr(level, isr, arg);
}
#endif

void LSM9DS0::uninstallISR(INTERRUPT_PINS_T intr)
{
  if (getPin(intr))
    {
      getPin(intr)->isrExit();
      delete getPin(intr);
      
      getPin(intr) = 0;
    }
}

mraa::Gpio*& LSM9DS0::getPin(INTERRUPT_PINS_T intr)
{
  switch(intr)
    {
    case INTERRUPT_G_INT:
      return m_gpioG_INT;
      break;
    case INTERRUPT_G_DRDY:
      return m_gpioG_DRDY;
      break;
    case INTERRUPT_XM_GEN1:
      return m_gpioXM_GEN1;
      break;
    case INTERRUPT_XM_GEN2:
      return m_gpioXM_GEN2;
      break;
    default:
      throw std::out_of_range(string(__FUNCTION__) +
                              ": Invalid interrupt enum passed");
    }
}
