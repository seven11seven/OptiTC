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
        self.f     = np.zeros((3,))    # rope
        self.norm_f = 1.0

        self.vtau_3 = 0.0

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

        self.btheta_3 = 0.0
        self.btheta_4 = 0.0
        self.btheta_5 = 0.0

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
    
    def compute_vtau_3(self):
        tmp_0 = self.h-self.theta_5+self.r
        tmp_1 = (self.m_4*self.btheta_5) / tmp_0
        tmp_2 = (self.m_4*self.atheta_5+self.m_4*self.g)*self.vtheta_5 / (tmp_0*tmp_0)
        self.vtau_3 = tmp_1 + tmp_2
        return self.vtau_3

    def compute_theta_2(self):
        tmp_1 = (self.m_4*self.atheta_4/self.tau_3+self.theta_4)^2
        tmp_2 = (self.m_4*self.atheta_3/self.tau_3+self.theta_3)^2
        self.theta_2 = np.power((tmp_1+tmp_2),0.5)
        return self.theta_2
    
    def compute_vtheta_2(self):
        tmp_1 = np.sqrt(self.theta_2)
        tmp_2 = self.m_4*self.atheta_4/self.tau_3 + self.theta_4
        tmp_3 = self.m_4*self.btheta_4/self.tau_3 - self.m_4*self.atheta_4*self.vtau_3/(self.tau_3*self.tau_3) + self.vtheta_4
        tmp_4 = self.m_4*self.atheta_3/self.tau_3 + self.theta_3
        tmp_5 = self.m_4*self.btheta_3/self.tau_3 - self.m_4*self.atheta_3*self.vtau_3/(self.tau_3*self.tau_3) + self.vtheta_3
        self.vtheta_2 = tmp_1*(tmp_2*tmp_3 + tmp_4*tmp_5)
        return self.vtheta_2

    def compute_theta_1(self):
        tmp_1 = (self.m_4*self.atheta_4)/(self.tau_3*self.theta_2) + self.theta_4/self.theta_2
        self.theta_1 = np.arccos(tmp_1)
        return self.theta_1

    def compute_vtheta_1(self):
        tmp_1 = (self.m_4*self.atheta_4)/(self.tau_3*self.theta_2) + self.theta_4/self.theta_2
        tmp_2 = -1 / (np.sqrt(1-tmp_1*tmp_1))
        tmp_3 = self.m_4*self.btheta_4/(self.tau_3*self.theta_2)
        tmp_4 = self.m_4*self.atheta_4*(self.vtau_3*self.theta_2+self.tau_3*self.vtheta_2) / np.square(self.tau_3*self.theta_2)
        tmp_5 = self.vtheta_4 / self.theta_2 - (self.theta_4*self.vtheta_2) / (self.theta_2*self.theta_2)
        self.vtheta_1 = tmp_2 * (tmp_3-tmp_4+tmp_5)
        return self.vtheta_1

    def compute_tau_1(self):
        tmp_1 = self.m_2*self.l*self.l*self.atheta_1 + self.Izz_2*self.atheta_1 + self.Ixx_3*self.atheta_1 + \
                self.m_3*self.atheta_1*self.theta_2*self.theta_2 + 2*self.m_3*self.theta_2*self.vtheta_1*self.vtheta_2
        tmp_2 = -self.theta_2*np.cos(self.theta_1) * (self.theta_3 + self.theta_2*np.sin(self.theta_1))
        tmp_3 = -self.theta_2*np.sin(self.theta_1) * (self.theta_4 - self.theta_2*np.cos(self.theta_1))
        self.tau_1 = tmp_1 + self.tau_3*(tmp_2+tmp_3)
        return self.tau_1
    
    def compute_tau_2(self):
        tmp_1 = self.m_3*self.atheta_2 - self.m_3*self.theta_2*self.vtheta_1*self.vtheta_1
        tmp_2 = -np.sin(self.theta_1) * (self.theta_3 + self.theta_2*np.sin(self.theta_1))
        tmp_3 = np.cos(self.theta_1) * (self.theta_4 - self.theta_2*np.cos(self.theta_1))
        self.tau_2 = tmp_1 + self.tau_3*(tmp_2+tmp_3)
        return self.tau_2
    
    def compute_f(self):
        tmp_1 = self.theta_3 + self.theta_2*np.sin(self.theta_1)
        tmp_2 = self.theta_4 - self.theta_2*np.cos(self.theta_1)
        tmp_3 = self.theta_5 - self.h + self.r
        self.f = self.tau_3 * np.ndarray([tmp_1, tmp_2, tmp_3])
        return self.f

    def compute_norm_f(self):
        self.norm_f = np.sqrt(np.sum(np.square(self.f)))
        return self.norm_f

    # How to compute atheata_1 and atheata_2 ???
    def compute_atau_3(self):
        """ Numerical Methods """
        pass

    def compute_atheta_2(self):
        """ Numerical Methods """
        pass

    def compute_atheta_1(self):
        """ Numerical Methods """
        pass

    def update(self):
        self.compute_tau_3()
        self.compute_vtau_3()
        self.compute_theta_2()
        self.compute_vtheta_2()
        self.compute_theta_1()
        self.compute_vtheta_1()
        
        # compute atheta_2
        # compute atheta_1
        self.compute_tau_1()
        self.compute_tau_2()
        self.compute_f()
        self.compute_norm_f()


if __name__ == "__main__":
    print(np.sqrt(np.sum(np.square([1, 1, 1]))))
