from tsBNgen import *
class tsBNgen:
    def __init__(self,T,N,N_level,Mat,Node_Type,CPD,Parent,CPD2,Parent2,loopbacks,CPD3=None,Parent3=None,loopbacks2=None,custom_time=0):
        '''
        A class to generate time series according to arbitrary dynamic Bayesian network structure.

        Attributes
        -----------
        T : int
            Length of each time series.
        
        N : int
            Number of time series.

        N_level : list
            Number of levels for the discrete nodes. Ignore this for the continuous nodes.

        Mat : data-frame
            Adjacency matrix corresponding to the Bayesian network at initial time.

        Node_Type : list
            Identifying nodes as either discrete "D" or continuous "C".

        CPD : dict
            Probability distribution fof the nodes at initial time point.

        Parent : dict
            Parents of each node at initial time point.

        CPD2 : dic
            Probability distribution of the nodes after initial time.

        Parent2 : dict
            Parents of each node for time points after the initial time point.
        
         loopbacks : dict
            Determining the temporal conection between the nodes.

        CPD3 : dict
            Probability distribution for the nodes. Use this entry when BN_sample_gen_loopback() is called.
                It defaults to empty.

        Parent3 : dict
            Parents of each node after the initial time point.  It is default to empty. 
                Use this entry when BN_sample_gen_loopback() is called 

        loopbacks2 : dict
            Determining the temporal conection between nodes. It is default to empty.
                Use this entry when BN_sample_gen_loopback() is called
        
        custom_time: int
            Determines at which time point, the new BN is used. The default is 0, which means
            the program learns it automatically from the loopbacks entry. 
            
        Methods
        ------------
        BFS(Row)
        Perform Breadth-first search for the given node(row).

        zero_loc(List)
            Find the index of zero values in a list.

        Role_Assignment()
            Identify the root node.
        
        DAG_ordering()
            Find the topological ordering of the graph.

        Child(Row)
            Finds the children of the node specified by the row of the adjacency matrix.  
        
        Multinomial_Select(index1,index2.ii)
            Generate sample according to the Multinomial distribution.
        
        Roots_length()
            Identifying the number of root nodes. 
        
        int_to_str(List)
            Concatenates list elements to a string.
        
        Valid_BN(parent)
            Verify whether the parent-child relationships between nodes are valid.
        
        Initial_sample() 
            Generate samples for all the nodes at initial time (t=0).
        
        Gaussian_select(sindex1,index2,ii=0)
            Generate sample according to the Gaussian distribution.
        
        continous_cpd()
            Identify which CPD entry to sample from.
        
        BN_sample():
            Generate samples for all the nodes after the initial time.
        
        BN_data_gen()
            Use this function under the following conditions:
                custom_time variable is not specified and the value of the loopback for all the variables is at most 1

        BN_sample_loopback()
            Generate samples for all the nodes for time t=k.
        
        BN_sample_gen_loopback()
            custom_time is not specified and you want the loopback value for some nodes to be at most 2.
                custom_time is specified and it is at least equal to the maximum loopback value of the loopbacks2.
        
        '''
        self.T=T
        self.N=N
        self.Mat=Mat
        self.Node_Type=Node_Type
        self.CPD=CPD 
        self.Node=[[] for ii in range(self.Mat.shape[0])]
        self.Parent=Parent 
        self.N_level=N_level
        self.CPD2=CPD2
        self.Parent2=Parent2
        self.flag=0
        if CPD3 is None:
            CPD3={}  
        self.CPD3=CPD3
        if Parent3 is None:
            Parent3={}
        self.Parent3=Parent3
        self.loopbacks=loopbacks 
        if loopbacks2 is None:
            loopbacks2={}
        self.loopbacks2=loopbacks2 
        self.custom_time=custom_time

    def BFS(self,Row):
        '''
        Perform Breadth-first search for the given node(row).

        Parameters
        --------------

        Row : int
              Corresponds to the row (node) in adjacency matrix.

        Returns
        --------------
        list
            The node and all its children.
        '''
        child = [Row]
        queue = []
        queue.extend(np.nonzero(self.Mat.iloc[Row, :].values)[0].tolist())
        while queue:
            vertex = queue.pop(0)
            if vertex not in child:
                child.append(vertex)
                queue.extend(np.nonzero(self.Mat.iloc[vertex, :].values)[0].tolist())
        return child

    @staticmethod  
    def zero_loc(List):
        '''
        Find the index of zero values in a list.

        Parameters
        ------------
        List: list

        Returns
        -----------
        list
            indices of zero values in a list.
        '''
        ind=[count  for count,ii in enumerate(List) if ii==0]
        return ind
    
    def __repr__(self):

        return  f''' length of each time series is {self.T}
        Number of time series samples is {self.N}
        The adjacency matrix is {self.Mat}
        Node Types are {self.Node_Type}
        Conditional Probability Table for initial time is {self.CPD}
        Parents of each node at initial time are {self.Parent}
        The number of levels for each discrete variable is {self.N_level}
        The roles are {self.Role}
        Conditional Probability Table for t2 to t_loopback is {self.CPD2}
        BN parents for time t2 ...t_loopback  for each node are {self.Parent2}
        Conditional Probability Table for t_loopback  to tn is {self.CPD3}
        BN parents for time t_loopback  ... tn for each node are {self.Parent3}
        '''

    def Role_Assignment(self):
        '''
        Identify the root node.

        Parameters
        ------------
        None

        Returns
        ------------
        None

        '''
        self.Role=[0]*len(self.top_order)
        for count,ii in enumerate(self.top_order):
            if(sum(self.Mat.iloc[:,ii]) ==0):
                self.Role[ii]=1  

    def DAG_ordering(self):
        '''
        Find the topological ordering of the graph

        Parameters
        -------------
        None

        Returns
        ------------
        None
        '''
        queue = []
        in_degree = np.count_nonzero(self.Mat, axis=0).tolist()
        index = tsBNgen.zero_loc(in_degree)  
        queue.extend(index)  
        visited_count = 0
        self.top_order = []
        while queue:
            P = queue.pop(0)
            Neighbor = self.Child(P)
            self.top_order.append(P)
            for count, ii in enumerate(Neighbor):
                in_degree[ii] = in_degree[ii] - 1
                if (in_degree[ii] == 0):
                    queue.append(ii)

            visited_count = visited_count + 1
        if (visited_count != self.Mat.shape[0]):
            print('DAG has a cycle')
    
    def Child(self, Row):  
        '''
        Finds the children of the node specified by the row of the adjacency matrix.

        Parameters
        -----------
        Row : int
              The row in the adjacency matrix, corresponding to the same node in a Bayesian network.

        Returns
        ----------
        list
             All the children of the given node.
        '''
        child=[]
        child.extend(np.nonzero(self.Mat.iloc[Row, :].values)[0].tolist())
        return child
    
    
    def Multinomial_Select(self,index1,index2,ii=0):
        '''
        Generate sample according to the Multinomial distribution.

        Parameters
        -----------
        index1: string
                key values of dictionary in CPD/CPD2/CPD3
        index2: int
                Determine which CPD entry to select.
        ii : int
             The node to generate the sample for. It defaults to 0.

        Returns
        ------------
        int
            The new generated sample.
        '''
        if(self.flag==0):
            if(self.Role[ii]==1):  
                Num=np.random.multinomial(1, self.CPD[str(index1)], size=1)
                Pos=Num.tolist()[0].index(1)
                return Pos+1
            else:  
                Num=np.random.multinomial(1, self.CPD[str(index1)][index2], size=1)
                Pos=Num.tolist()[0].index(1)
                return Pos+1
        elif (self.flag==1):
            if(len(self.Parent2[str(ii)])!=0):
                Num=np.random.multinomial(1, self.CPD2[str(index1)][index2], size=1)
                Pos=Num.tolist()[0].index(1)
            else:
                Num=np.random.multinomial(1, self.CPD2[str(index1)], size=1)
            Pos=Num.tolist()[0].index(1)    
            return Pos+1
        
        elif (self.flag==2):
            if(len(self.Parent3[str(ii)])!=0):
                Num=np.random.multinomial(1, self.CPD3[str(index1)][index2], size=1)
            else:
                Num=np.random.multinomial(1, self.CPD3[str(index1)], size=1)
            Pos=Num.tolist()[0].index(1)
            return Pos+1

    
    
    def parents_len(self,Node): 
        return np.count_nonzero(self.Mat.iloc[:,Node], axis=0).tolist() 
    
    def Roots_length(self): 
        '''
        Identifying the number of root nodes. 

        Parameters
        ------------
        None

        Returns
        ------------
        int
            Number of root nodes.
        '''
        parent=np.count_nonzero(self.Mat, axis=0).tolist()  
        return len([idx  for idx,ii in enumerate(parent) if ii==0])

    @staticmethod
    def int_to_str(List):
        '''
        Concatenates list elements to a string.

        Parameters
        -------------
        List : list

        Returns
        -------------
        string
              concatenated list elements as a string.

        '''
        if not isinstance(List, list): 
            List=list(map(int, str(List)))
        List=[str(ii) for ii in List]
        return ''.join(List)
    
    def Valid_BN(self,parent): 
        '''
        Verify whether the parent-child relationships between nodes are valid.

        Parameters
        ------------
        parent : dict
            dictionary where the keys are the nodes and the values are the list of parents.

        Returns
        -----------
        None

        Raises
        -----------
        Exception
            "Parent of a discrete node cannot be continuous"
        '''
        for count,ii in enumerate(self.top_order):
            if(self.Node_Type[ii]=='D' and any(self.Node_Type[jj]=='C' for jj in parent[str(ii)])):
                raise Exception("Parent of a discrete node cannot be continuous")

    
    def Initial_sample(self): 
        '''
        Generate samples for all the nodes at initial time (t=0)

        Parameters
        --------------
        None

        Returns
        -------------
        None

        Raises
        ------------
        Exception
            Parent of a discrete node cannot be continuous
        '''
        self.DAG_ordering()
        self.Role_Assignment()
        self.flag=0
        self.Valid_BN(self.Parent)

        for count,ii in enumerate(self.top_order):
            parent=tsBNgen.int_to_str(self.Parent[str(ii)])+str(ii)  
            if(self.Role[ii]==1 and self.Node_Type[ii]=='D'):
                self.Node[ii].append(self.Multinomial_Select(parent,0,ii))   
            elif(self.Role[ii]==1 and self.Node_Type[ii]=='C'):
                self.Node[ii].extend(self.Gaussian_select(ii,0)) 
            elif(all(self.Node_Type[jj]=='D' for jj in self.Parent[str(ii)])): 
                all_parents=[]
                parent_N_level=[]
                for count2,jj in enumerate(self.Parent[str(ii)]):
                    all_parents.append(self.Node[jj][-1])
                    parent_N_level.append(self.N_level[jj])
           
                self.all_parents=all_parents  
                self.parent_N_level=parent_N_level 
                CPD_entry=self.continous_cpd()

                if(self.Node_Type[ii]=='D'):
                    self.Node[ii].append(self.Multinomial_Select(parent,CPD_entry,ii))
                elif(self.Node_Type[ii]=='C'):
                    self.Node[ii].extend(self.Gaussian_select(parent,CPD_entry))

            elif(any(self.Node_Type[jj]=='D' for jj in self.Parent[str(ii)])):
                D_Parent=[kk for kk in self.Parent[str(ii)] if self.Node_Type[kk]=='D'] 
                all_parents=[]
                parent_N_level=[]
                for count2,jj in enumerate(D_Parent):
                    all_parents.append(self.Node[jj][-1])
                    parent_N_level.append(self.N_level[jj])
                
                self.all_parents=all_parents
                self.parent_N_level=parent_N_level
                CPD_entry=self.continous_cpd()
                
                C_Parent=[kk for kk in self.Parent[str(ii)] if self.Node_Type[kk]=='C'] 
                temp=0
                for count3,kk in enumerate(C_Parent):  
                    cont_SUM=self.CPD[str(parent)][str(kk)+str(ii)]['coefficient'][0][CPD_entry]*self.Node[kk][-1]
                    temp=temp+cont_SUM
                intercept=np.random.normal(0,self.CPD[str(parent)]['sigma_intercept'][CPD_entry],1)
                self.Node[ii].extend(self.Gaussian_select(temp+intercept,self.CPD[str(parent)]['sigma'][CPD_entry],ii))

            elif(all(self.Node_Type[jj]=='C' for jj in self.Parent[str(ii)])):  
                temp=0
                for count3,kk in enumerate(self.Parent[str(ii)]):  
                    cont_SUM=self.CPD[str(parent)][str(kk)+str(ii)]['coefficient'][0][0]*self.Node[kk][-1]
                    temp=temp+cont_SUM
                intercept=np.random.normal(0,self.CPD[str(parent)]['sigma_intercept'][0],1)
                self.Node[ii].extend(self.Gaussian_select(temp+intercept,self.CPD[str(parent)]['sigma'][0],ii))

    
    def Gaussian_select(self,index1,index2,ii=0):   
        if(self.flag==0):
            C_Parent=[kk for kk in self.Parent[str(ii)] if self.Node_Type[kk]=='C']
            if(len(C_Parent)!=0):
                return (np.random.normal(index1,index2, 1).tolist())
            else:
                return (np.random.normal(self.CPD[str(index1)]['mu'+str(index2)],self.CPD[str(index1)]['sigma'+str(index2)], 1).tolist())
        
        elif(self.flag==1):
            C_Parent=[kk for kk in self.Parent2[str(ii)] if self.Node_Type[kk]=='C']
            if(len(C_Parent)!=0):
                return (np.random.normal(index1,index2, 1).tolist())
            else:
                return (np.random.normal(self.CPD2[str(index1)]['mu'+str(index2)],self.CPD2[str(index1)]['sigma'+str(index2)], 1).tolist())
        elif(self.flag==2):
            C_Parent=[kk for kk in self.Parent3[str(ii)] if self.Node_Type[kk]=='C']
            if(len(C_Parent)!=0):
                return (np.random.normal(index1,index2, 1).tolist())
            else:
                return (np.random.normal(self.CPD3[str(index1)]['mu'+str(index2)],self.CPD3[str(index1)]['sigma'+str(index2)], 1).tolist())

    
    def Level_multiplied(self):
        self.level_multiply=[]
        Temp=self.parent_N_level  
        for ii in range(len(Temp)):
            A1=reduce(lambda x, y: x*y,Temp)
            self.level_multiply.append(A1)
            Temp=Temp[1:]
        self.level_multiply=self.level_multiply[1:]
        self.level_multiply.insert(len(self.parent_N_level ),1)

    def continous_cpd(self):
        self.Level_multiplied()
        List=[(ii-1)*jj for ii, jj in zip(self.all_parents,self.level_multiply)]
        return sum(List)

    def BN_sample(self):
        '''
        Generate samples for all the nodes after the initial time.

        Parameters
        --------------
        None

        Returns
        -------------
        None

        Notes
        ------------
        Use this function to generate samples if the loopback values are at most one. 
        Loopback=1 means that a node at time t is connected to the node at t-1.

        Raises
        ------------
        Exception
            Parent of a discrete node cannot be continuous
        '''
        self.flag=1
        loopbacks_temp=self.loopbacks.copy()
        self.Valid_BN(self.Parent2)
        for count,ii in enumerate(self.top_order):
            parent=tsBNgen.int_to_str(self.Parent2[str(ii)])+str(ii) 
            if(all(self.Node_Type[jj]=='D' for jj in self.Parent2[str(ii)])):
                temp=0
                all_parents=[]
                parent_N_level=[] 
                for count2,jj in enumerate(self.Parent2[str(ii)]): 
                    loopbacks_keys=tsBNgen.int_to_str(jj)+str(ii)  
                    if(loopbacks_keys in self.loopbacks.keys()): 
                        if(self.top_order.index(jj) < count):  
                            for count3,mm in enumerate(self.loopbacks[loopbacks_keys]): 
                                all_parents.append(self.Node[jj][-1-mm]) 
                        else: 
                            for count3,mm in enumerate(self.loopbacks[loopbacks_keys]): 
                                LP=mm-1 
                                all_parents.append(self.Node[jj][-1-LP]) 
                        del self.loopbacks[loopbacks_keys]
                    elif(jj != ii): 
                        all_parents.append(self.Node[jj][-1])  
                    parent_N_level.append(self.N_level[jj])

                self.loopbacks=loopbacks_temp.copy()
               
                self.all_parents=all_parents  
                self.parent_N_level=parent_N_level 
                CPD_entry=self.continous_cpd()
                
                if(self.Node_Type[ii]=='D'):
                    self.Node[ii].append(self.Multinomial_Select(parent,CPD_entry,ii))
                elif(self.Node_Type[ii]=='C'):
                    self.Node[ii].extend(self.Gaussian_select(parent,CPD_entry))
            
            elif(any(self.Node_Type[jj]=='D' for jj in self.Parent2[str(ii)])):
                loopbacks_temp=self.loopbacks.copy()
                D_Parent=[kk for kk in self.Parent2[str(ii)] if self.Node_Type[kk]=='D'] 
                all_parents=[]
                parent_N_level=[]
                for count2,jj in enumerate(D_Parent):
                    loopbacks_keys=tsBNgen.int_to_str(jj)+str(ii)  
                    if(loopbacks_keys in self.loopbacks.keys()): 
                        if(self.top_order.index(jj) < count):  
                            for count3,mm in enumerate(self.loopbacks[loopbacks_keys]): 
                                all_parents.append(self.Node[jj][-1-mm]) 
                        else: 
                            for count3,mm in enumerate(self.loopbacks[loopbacks_keys]): 
                                LP=mm-1 
                                all_parents.append(self.Node[jj][-1-LP]) 
                        del self.loopbacks[loopbacks_keys]
                    elif(jj != ii): 

                        all_parents.append(self.Node[jj][-1])  
                    parent_N_level.append(self.N_level[jj])

                    self.loopbacks=loopbacks_temp.copy()
                self.all_parents=all_parents
                self.parent_N_level=parent_N_level
                CPD_entry=self.continous_cpd()
                
                C_Parent=[kk for kk in self.Parent2[str(ii)] if self.Node_Type[kk]=='C'] 
                temp=0
                loopbacks_temp=self.loopbacks.copy()
                for count3,kk in enumerate(C_Parent):
                   loopbacks_keys=tsBNgen.int_to_str(kk)+str(ii)
                   if(loopbacks_keys in self.loopbacks.keys()): 
                       if(self.top_order.index(kk) < count): 
                           for count3,mm in enumerate(self.loopbacks[loopbacks_keys]):
                              cont_SUM=self.CPD2[str(parent)][str(kk)+str(ii)]['coefficient'][mm][CPD_entry]*self.Node[kk][-1-mm]
                              temp=temp+cont_SUM
                       else: 
                           for count3,mm in enumerate(self.loopbacks[loopbacks_keys]): 
                               LP=mm-1 
                               cont_SUM=self.CPD2[str(parent)][str(kk)+str(ii)]['coefficient'][LP][CPD_entry]*self.Node[kk][-1-LP]
                               temp=temp+cont_SUM
                       del self.loopbacks[loopbacks_keys]
                   elif(kk != ii):
                       cont_SUM=self.CPD2[str(parent)][str(kk)+str(ii)]['coefficient'][0][CPD_entry]*self.Node[kk][-1]
                       temp=temp+cont_SUM
                intercept=np.random.normal(0,self.CPD2[str(parent)]['sigma_intercept'][CPD_entry],1)
                self.Node[ii].extend(self.Gaussian_select(temp+intercept,self.CPD2[str(parent)]['sigma'][CPD_entry],ii))
                self.loopbacks=loopbacks_temp.copy()
            
            elif(all(self.Node_Type[jj]=='C' for jj in self.Parent2[str(ii)])):
                C_Parent=[kk for kk in self.Parent2[str(ii)] if self.Node_Type[kk]=='C'] 
                temp=0
                loopbacks_temp=self.loopbacks.copy()
                for count3,kk in enumerate(C_Parent):
                   loopbacks_keys=tsBNgen.int_to_str(kk)+str(ii)
                   if(loopbacks_keys in self.loopbacks.keys()): 
                       if(self.top_order.index(kk) < count):  
                           for count3,mm in enumerate(self.loopbacks[loopbacks_keys]):
                              cont_SUM=self.CPD2[str(parent)][str(kk)+str(ii)]['coefficient'][mm][0]*self.Node[kk][-1-mm]
                              temp=temp+cont_SUM
                       else: 
                           for count3,mm in enumerate(self.loopbacks[loopbacks_keys]): 
                               LP=mm-1 
                               cont_SUM=self.CPD2[str(parent)][str(kk)+str(ii)]['coefficient'][LP][0]*self.Node[kk][-1-LP]
                               temp=temp+cont_SUM
                       del self.loopbacks[loopbacks_keys]
                   elif(kk != ii):
                       cont_SUM=self.CPD2[str(parent)][str(kk)+str(ii)]['coefficient'][0][0]*self.Node[kk][-1]
                       temp=temp+cont_SUM
                intercept=np.random.normal(0,self.CPD2[str(parent)]['sigma_intercept'][0],1)
                self.Node[ii].extend(self.Gaussian_select(temp+sintercept,self.CPD2[str(parent)]['sigma'][0],ii))
                self.loopbacks=loopbacks_temp.copy()
            

        
    def BN_data_gen(self):
        '''
        It uses Initial_sample for initial time(t=0) and BN_sample for time point t=1 up to time t=T (length of time series)

        Parameters
        -------------
        None

        Returns
        -------------
        None
            None

        Raises
        ------------
        Exception
            Parent of a discrete node cannot be continuous

        Notes
        -------------
        Use this function under the following conditions: custom_time variable is not specified 
        and the value of the loopback for all the variables is at most 1
        '''
        keys = range(len(self.Node)) 
        self.BN_Nodes= dict(zip(keys, ([[] for ii in range(self.N)] for _ in keys )))
        for ii in range(self.N):
            
            self.Initial_sample()
          
            for kk in range(1,self.T):
                self.BN_sample()
            for jj in range(len(self.Node)):
                self.BN_Nodes[jj][ii]=self.Node[jj] 
                self.Node[jj]=[]                    
              
        
    
    def BN_sample_loopback(self):
        '''
        Generate samples for all the nodes given CPD3 and Parent3 are used.

        Parameters
        --------------
        None

        Returns
        -------------
        None

        Raises
        ------------
        Exception
            Parent of a discrete node cannot be continuous
        
        Notes
        ------------
        Use this function when you want to incorporate three BNs.
        '''
        self.flag=2
        self.Valid_BN(self.Parent3)
        self.DAG_ordering()  
        loopbacks_temp=self.loopbacks2.copy()  
        for count,ii in enumerate(self.top_order):  
            parent=tsBNgen.int_to_str(self.Parent3[str(ii)])+str(ii) 
            if(all(self.Node_Type[jj]=='D' for jj in self.Parent3[str(ii)])):
                temp=0
                all_parents=[]
                parent_N_level=[]  
                for count2,jj in enumerate(self.Parent3[str(ii)]): 
                    loopbacks_keys=tsBNgen.int_to_str(jj)+str(ii)  
                    if(loopbacks_keys in self.loopbacks2.keys()): 
                        if(self.top_order.index(jj) < count):  
                            for count3,mm in enumerate(self.loopbacks2[loopbacks_keys]): 
                                all_parents.append(self.Node[jj][-1-mm]) 
                        else: 
                            for count3,mm in enumerate(self.loopbacks2[loopbacks_keys]): 
                                LP=mm-1 
                                all_parents.append(self.Node[jj][-1-LP]) 
                        del self.loopbacks2[loopbacks_keys]

                    
                    elif(jj != ii): 
                        all_parents.append(self.Node[jj][-1])  
                    parent_N_level.append(self.N_level[jj])

                self.loopbacks2=loopbacks_temp.copy()
                self.all_parents=all_parents  
                self.parent_N_level=parent_N_level  
                CPD_entry=self.continous_cpd()
                if(self.Node_Type[ii]=='D'):
                    self.Node[ii].append(self.Multinomial_Select(parent,CPD_entry,ii))
                elif(self.Node_Type[ii]=='C'):
                    self.Node[ii].extend(self.Gaussian_select(parent,CPD_entry))
            elif(any(self.Node_Type[jj]=='D' for jj in self.Parent3[str(ii)])):
                loopbacks_temp=self.loopbacks2.copy()
                D_Parent=[kk for kk in self.Parent3[str(ii)] if self.Node_Type[kk]=='D']
                temp=0
                all_parents=[]
                parent_N_level=[]  
                for count2,jj in enumerate(D_Parent):  
                    loopbacks_keys=tsBNgen.int_to_str(jj)+str(ii)  
                    if(loopbacks_keys in self.loopbacks2.keys()): 
                        if(self.top_order.index(jj) < count):  
                            for count3,mm in enumerate(self.loopbacks2[loopbacks_keys]): 
                                all_parents.append(self.Node[jj][-1-mm]) 
                        else: 
                            for count3,mm in enumerate(self.loopbacks2[loopbacks_keys]): 
                                LP=mm-1 
                                all_parents.append(self.Node[jj][-1-LP]) 
                        del self.loopbacks2[loopbacks_keys]

                    
                    elif(jj != ii): 
                        all_parents.append(self.Node[jj][-1])  
                    parent_N_level.append(self.N_level[jj])

                self.loopbacks2=loopbacks_temp.copy()
                self.all_parents=all_parents  
                self.parent_N_level=parent_N_level  
                CPD_entry=self.continous_cpd()

                C_Parent=[kk for kk in self.Parent3[str(ii)] if self.Node_Type[kk]=='C'] 
                temp=0
                loopbacks_temp=self.loopbacks2.copy()
                for count3,kk in enumerate(C_Parent):
                   loopbacks_keys=tsBNgen.int_to_str(kk)+str(ii)
                   if(loopbacks_keys in self.loopbacks2.keys()): 
                       if(self.top_order.index(kk) < count):  
                           for count3,mm in enumerate(self.loopbacks2[loopbacks_keys]):
                              cont_SUM=self.CPD3[str(parent)][str(kk)+str(ii)]['coefficient'][mm][CPD_entry]*self.Node[kk][-1-mm]
                              temp=temp+cont_SUM
                       else: 
                           for count3,mm in enumerate(self.loopbacks2[loopbacks_keys]): 
                               LP=mm-1 
                               cont_SUM=self.CPD3[str(parent)][str(kk)+str(ii)]['coefficient'][LP][CPD_entry]*self.Node[kk][-1-LP]
                               temp=temp+cont_SUM
                       del self.loopbacks[loopbacks_keys]
                   elif(kk != ii):
                       cont_SUM=self.CPD3[str(parent)][str(kk)+str(ii)]['coefficient'][0][CPD_entry]*self.Node[kk][-1]
                       temp=temp+cont_SUM
                intercept=np.random.normal(0,self.CPD3[str(parent)]['sigma_intercept'][CPD_entry],1)
                self.Node[ii].extend(self.Gaussian_select(temp+intercept,self.CPD3[str(parent)]['sigma'][CPD_entry],ii))
                self.loopbacks2=loopbacks_temp.copy()
            
            elif(all(self.Node_Type[jj]=='C' for jj in self.Parent3[str(ii)])):
                C_Parent=[kk for kk in self.Parent3[str(ii)] if self.Node_Type[kk]=='C'] 
                temp=0
                loopbacks_temp=self.loopbacks2.copy()
                for count3,kk in enumerate(C_Parent):
                   loopbacks_keys=tsBNgen.int_to_str(kk)+str(ii)
                   if(loopbacks_keys in self.loopbacks2.keys()): 
                       if(self.top_order.index(kk) < count):  
                           for count3,mm in enumerate(self.loopbacks2[loopbacks_keys]):
                              cont_SUM=self.CPD3[str(parent)][str(kk)+str(ii)]['coefficient'][mm][0]*self.Node[kk][-1-mm]
                              temp=temp+cont_SUM
                       else: 
                           for count3,mm in enumerate(self.loopbacks2[loopbacks_keys]): 
                               LP=mm-1 
                               cont_SUM=self.CPD3[str(parent)][str(kk)+str(ii)]['coefficient'][LP][0]*self.Node[kk][-1-LP]
                               temp=temp+cont_SUM
                       del self.loopbacks[loopbacks_keys]
                   elif(kk != ii):
                       cont_SUM=self.CPD3[str(parent)][str(kk)+str(ii)]['coefficient'][0][0]*self.Node[kk][-1]
                       temp=temp+cont_SUM
                intercept=np.random.normal(0,self.CPD3[str(parent)]['sigma_intercept'][0],1)
                self.Node[ii].extend(self.Gaussian_select(temp+intercept,self.CPD3[str(parent)]['sigma'][0],ii))
                self.loopbacks2=loopbacks_temp.copy()

    
    def BN_sample_gen_loopback(self):
        '''
        Generate time series data for all the nodes for all time. See Notes for the more information.

        Parameters
        -------------
        None

        Returns
        -------------
        None

        Raises
        ------------
        Exception
            Parent of a discrete node cannot be continuous
        
        Notes
        ------------
        This is more general form of BN_data_gen that supports only two different BN structures 
        or loopback value of maximum one for all the nodes.
        '''
        keys = range(len(self.Node)) 
        self.BN_Nodes= dict(zip(keys, ([[] for ii in range(self.N)] for _ in keys )))
        Max_loopback=max(sum(self.loopbacks2.values(),[]))

        if (self.custom_time == 0):
            for ii in range(self.N):
                self.Initial_sample()
                for kk in range(1,Max_loopback):
                    self.BN_sample()
                for mm in range(Max_loopback,self.T):
                    self.BN_sample_loopback()
                for jj in range(len(self.Node)):
                    self.BN_Nodes[jj][ii]=self.Node[jj] 
                    self.Node[jj]=[]    
        else:
            for ii in range(self.N):
                self.Initial_sample()
                for kk in range(1,self.custom_time):
                    self.BN_sample()
                for mm in range(self.custom_time,self.T):
                    self.BN_sample_loopback()
                for jj in range(len(self.Node)):
                    self.BN_Nodes[jj][ii]=self.Node[jj] 
                    self.Node[jj]=[]

              



      

         




