import numpy as np

class Config:
    """
    Configuration for the Crane-Module System
    """
    def __init__(self) -> None:
        ## Temporal variables
        self.tau_1 = 1.0    # rotation joint
        self.tau_2 = 1.0    # linear joint
        self.tau_3 = 1.0    # f/l
        self.f     = 1.0    # rope

        self.theta_1 = 1.0  # rotation joint
        self.theta_2 = 1.0  # linear joint
        self.theta_3 = 1.0  # x coordinate of module
        self.theta_4 = 1.0  # y coordinate of module
        self.theta_5 = 1.0  # z coordinate of module

        self.vtheta_1 = 1.0  # rotation joint
        self.vtheta_2 = 1.0  # linear joint
        self.vtheta_3 = 1.0  # dx/dt
        self.vtheta_4 = 1.0  # dy/dt
        self.vtheta_5 = 1.0  # dz/dt

        self.atheta_1 = 1.0  # rotation joint
        self.atheta_2 = 1.0  # linear joint
        self.atheta_3 = 1.0  # ddx/dt^2
        self.atheta_4 = 1.0  # ddy/dt^2
        self.atheta_5 = 1.0  # ddz/dt^2

        ## Constant Parameters
        self.h = 1.0    # hight of the jib
        self.r = 1.0    # from module's CoM to hook point
        self.l = 1.0    # from hook point to trolley
        self.g = 9.8    # absolute value of gravity

        self.Izz_2 = 1.0  # inertia of jib
        self.Ixx_3 = 1.0  # inertia of trolley
        
        self.m_1 = 1.0  # mass of mast
        self.m_2 = 1.0  # mass of jib
        self.m_3 = 1.0  # mass of trolley
        self.m_4 = 1.0  # mass of module

    def compute_tau_3(self):
        self.tau_3 = (self.m_4*self.atheta_5+self.m_4*self.g)/(self.h-self.theta_5+self.r)
        return self.tau_3
    
    def compute_theta_1(self):
        # self.theta_2 = self.compute_theta_2()
        tmp_1 = (self.m_4*self.atheta_4)/(self.tau_3*self.theta_2) + self.theta_4/self.theta_2
        self.theta_1 = np.arccos(tmp_1)
        return self.theta_1

    def compute_theta_2(self):
        # self.tau_3 = self.compute_tau_3()
        tmp_1 = (self.m_4*self.atheta_4/self.tau_3+self.theta_4)^2
        tmp_2 = (self.m_4*self.atheta_3/self.tau_3+self.theta_3)^2
        self.theta_2 = np.power((tmp_1+tmp_2),0.5)
        return self.theta_2

    def compute_tau_1(self):
        pass