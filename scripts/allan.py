#!/usr/bin/env python
import rospy
import sys
import allantools
import rosbag
import numpy as np
import csv
import rospkg
import os
import matplotlib.pyplot as plt  # only for plotting, not required for calculations
import math

def getWhiteNoiseSegment(tau,sigma):

    m = -0.5 # slope of random walk
    """""""""""""""""""""""""""""""""
    " Find point where slope = -0.5 "
    """""""""""""""""""""""""""""""""

    i = 1
    idx = 1
    mindiff = 999
    logTau = -999
    while (logTau<0):
        logTau = math.log(tau[i],10)
        logSigma = math.log(sigma[i],10)
        prevLogTau = math.log(tau[i-1],10)
        prevLogSigma = math.log(sigma[i-1],10)
        slope = (logSigma-prevLogSigma)/(logTau-prevLogTau)
        diff = abs(slope-m)
        if (diff<mindiff):
            mindiff = diff
            idx = i
        i = i + 1

    """"""""""""""""""""""""""""""
    " Project line to tau = 10^0 "
    """"""""""""""""""""""""""""""
    x1 = math.log(tau[idx],10)
    y1 = math.log(sigma[idx],10)
    x2 = 0
    y2 = m*(x2-x1)+y1

    return (pow(10,x1),pow(10,y1),pow(10,x2),pow(10,y2))


def GetWhiteNoise(tau,sigma):

    m = -0.5 # slope of random walk
    sigma_min = getBiasInstabilityPoint(tau,sigma)
    window_size = 100
    log_sigma_avg = 0.0

    while(sigma_min[2]-window_size<60):
        window_size = window_size - 20

    for i in range(sigma_min[2]-window_size):
        log_tau = math.log(tau[i],10)
        log_sigma = math.log(sigma[i],10)
        log_sigma_avg = log_sigma_avg + log_sigma - m*log_tau
    log_sigma_avg = log_sigma_avg/float(sigma_min[2]-window_size)

    return (tau[0],pow(10,m*math.log(tau[0],10)+log_sigma_avg),
            tau[sigma_min[2]],pow(10,m*math.log(tau[sigma_min[2]],10)+log_sigma_avg),
            1,pow(10,log_sigma_avg))

def GetRandomWalk(tau,sigma):

    m = 0.5 # slope of random walk
    sigma_min = getBiasInstabilityPoint(tau,sigma)
    window_size = 100
    log_sigma_avg = 0.0

    while(sigma_min[2]+window_size>tau.shape[0]-60):
        window_size = window_size - 20

    for i in range(sigma_min[2]+window_size,tau.shape[0]):
        log_tau = math.log(tau[i],10)
        log_sigma = math.log(sigma[i],10)
        log_sigma_avg = log_sigma_avg + log_sigma - m*log_tau
    log_sigma_avg = log_sigma_avg/float(-sigma_min[2]-window_size+tau.shape[0])

    return (3,pow(10,m*math.log(3,10)+log_sigma_avg),
            tau[tau.shape[0]-1],pow(10,m*math.log(tau[tau.shape[0]-1],10)+log_sigma_avg),
            3,pow(10,m*math.log(3,10)+log_sigma_avg))


def getBiasInstabilityPoint(tau,sigma):
    i = 1
    while (i<tau.size):
        if (tau[i]>1) and ((sigma[i]-sigma[i-1])>0): # only check for tau > 10^0
            break
        i = i + 1
    return (tau[i],sigma[i],i)

def main(args):

    rospy.init_node('allan_variance_node')

    t0 = rospy.get_time()

    """"""""""""""
    " Parameters "
    """"""""""""""
    bagfile = rospy.get_param('~bagfile_path','~/data/static.bag')
    topic = rospy.get_param('~imu_topic_name','/imu')
    axis = rospy.get_param('~axis',0)
    sampleRate = rospy.get_param('~sample_rate',100)
    isDeltaType = rospy.get_param('~delta_measurement',False)
    numTau = rospy.get_param('~number_of_lags',1000)
    resultsPath = rospy.get_param('~results_directory_path',None)

    """"""""""""""""""""""""""
    " Results Directory Path "
    """"""""""""""""""""""""""
    if resultsPath is None:
        paths = rospkg.get_ros_paths()
        path = paths[1] # path to workspace's devel
        idx = path.find("ws/")
        workspacePath = path[0:(idx+3)]
        resultsPath = workspacePath + 'av_results/'

        if not os.path.isdir(resultsPath):
            os.mkdir(resultsPath)

    print "\nResults will be save in the following directory: \n\n\t %s\n"%resultsPath

    """"""""""""""""""
    " Form Tau Array "
    """"""""""""""""""
    taus = [None]*numTau

    cnt = 0;
    for i in np.linspace(-2.0, 5.0, num=numTau): # lags will span from 10^-2 to 10^5, log spaced
        taus[cnt] = pow(10,i)
        cnt = cnt + 1

    """""""""""""""""
    " Parse Bagfile "
    """""""""""""""""
    bag = rosbag.Bag(bagfile)

    N = bag.get_message_count(topic) # number of measurement samples

    data = np.zeros( (6,N) ) # preallocate vector of measurements

    if isDeltaType:
        scale = sampleRate
    else:
        scale = 1.0

    cnt = 0
    for topic, msg, t in bag.read_messages(topics=[topic]):
        data[0,cnt] = msg.linear_acceleration.x * scale
        data[1,cnt] = msg.linear_acceleration.y * scale
        data[2,cnt] = msg.linear_acceleration.z * scale
        data[3,cnt] = msg.angular_velocity.x * scale
        data[4,cnt] = msg.angular_velocity.y * scale
        data[5,cnt] = msg.angular_velocity.z * scale
        cnt = cnt + 1

    bag.close()

    print "[%0.2f seconds] Bagfile parsed\n"%(rospy.get_time()-t0)

    """"""""""""""""""
    " Allan Variance "
    """"""""""""""""""
    if axis is 0:
        currentAxis = 1 # loop through all axes 1-6
    else:
        currentAxis = axis # just loop one time and break


    fsummary = open(resultsPath+"allan_deviation.yaml","wt")
    while (currentAxis <= 6):
        (taus_used, adev, adev_err, adev_n) = allantools.oadev(data[currentAxis-1], data_type='freq', rate=float(sampleRate), taus=np.array(taus) )

        #WhiteNoiseSegment = getWhiteNoiseSegment(taus_used,adev)
        WhiteNoiseSegment = GetWhiteNoise(taus_used,adev) # white noise slope
        biasInstabilityPoint = getBiasInstabilityPoint(taus_used,adev)
        RandomWalkSegment = GetRandomWalk(taus_used,adev)

        biasInstability = biasInstabilityPoint[1]
        
        """""""""""""""
        " Save as CSV "
        """""""""""""""
        if (currentAxis==1):
            fname = 'allan_accel_x'
            title = 'Allan Deviation: Accelerometer X'
            axis_name = 'acc_x'
        elif (currentAxis==2):
            fname = 'allan_accel_y'
            title = 'Allan Deviation: Accelerometer Y'
            axis_name = 'acc_y'
        elif (currentAxis==3):
            fname = 'allan_accel_z'
            title = 'Allan Deviation: Accelerometer Z'
            axis_name = 'acc_z'
        elif (currentAxis==4):
            fname = 'allan_gyro_x'
            title = 'Allan Deviation: Gyroscope X'
            axis_name = 'gyr_x'
        elif (currentAxis==5):
            fname = 'allan_gyro_y'
            title = 'Allan Deviation: Gyroscope Y'
            axis_name = 'gyr_y'
        elif (currentAxis==6):
            fname = 'allan_gyro_z'
            title = 'Allan Deviation: Gyroscope Z'
            axis_name = 'gyr_z'

        print "[%0.2f seconds] Finished calculating allan variance - writing results to %s"%(rospy.get_time()-t0,fname)

        #f = open(resultsPath + fname + '.csv', 'wt')


        try:
            fsummary.write(axis_name+"_white_noise: "+str(WhiteNoiseSegment[5])+"\n") #white noise
            fsummary.write(axis_name+"_random_walk: "+str(RandomWalkSegment[5])+"\n") #random walk
            """
            writer = csv.writer(f)
            writer.writerow( ('Random Walk', 'Bias Instability') )
            writer.writerow( (WhiteNoiseSegment[5], biasInstability) )
            writer.writerow( ('Tau', 'AllanDev', 'AllanDevError', 'AllanDevN') )
            for i in range(taus_used.size):
                writer.writerow( (taus_used[i],adev[i],adev_err[i],adev_n[i])  )
            """

        finally:
            pass
            #f.close()

        """""""""""""""
        " Plot Result "
        """""""""""""""
        plt.figure(figsize=(12,8))
        ax = plt.gca()
        ax.set_yscale('log')
        ax.set_xscale('log')

        plt.plot(taus_used,adev)
        plt.plot([WhiteNoiseSegment[0],WhiteNoiseSegment[2]],
                 [WhiteNoiseSegment[1],WhiteNoiseSegment[3]],'r--')
        plt.plot(WhiteNoiseSegment[4],WhiteNoiseSegment[5],'rx',markeredgewidth=2.5,markersize=14.0)

        plt.plot([RandomWalkSegment[0],RandomWalkSegment[2]],
                 [RandomWalkSegment[1],RandomWalkSegment[3]],'b--')
        plt.plot(RandomWalkSegment[4],RandomWalkSegment[5],'bx',markeredgewidth=2.5,markersize=14.0)

        plt.plot(biasInstabilityPoint[0],biasInstabilityPoint[1],'go')

        plt.grid(True, which="both")
        plt.title(title)
        plt.xlabel('Tau (s)')
        plt.ylabel('ADEV')

        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                    ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)

        #plt.show(block=False)

        plt.savefig(resultsPath + fname)

        currentAxis = currentAxis + 1 + axis*6 # increment currentAxis also break if axis is not =0

    fsummary.close()
    inp=raw_input("Press Enter key to close figures and end program\n")

if __name__ == '__main__':
  main(sys.argv)
