module Main where

import Text.Printf (printf)
import System.Environment (getArgs)
import System.Process (readProcess)
import System.CPUTime (getCPUTime)
import Data.List (intercalate)

--
-- Config
--

-- Each jobs submitted to the cluster will run a "batch" of parameters
-- one-at-a-time.  Hopefully this avoids a little bit of overhead from
-- cluster scheduling.  So far, I've aimed to make this target ~ 30 min
-- jobs; still good for the short queue not too short.
batch::Int
batch = 5

params = [(tm, tau, bg) | tm <-[200..400]
                        , tau<-[5,10,15]
                        , bg <-[50,100,150,200,250]]   

-- recommend making a link to the bin directory from the cwd 
-- This makes the output text a little simpler to visually parse and avoids "baked in" paths
exe = "bin/ProcessStack_woGPU" 

-- This setup function forms the input file name
source :: Int -> String
source tm = printf "/nobackup/keller/clackn/mmu/14-05-21/Mmu_E1_CAGTAG1.TM%06d_timeFused_blending/SPM00_TM%06d_CM00_CM01_CHN00.fusedStack.corrected" tm tm 

--
--
--

data ParameterSet = ParameterSet { paramId :: Int
                                 , paramTM :: Int
                                 , paramTau :: Int
                                 , paramBackground :: Int
                                 }
  deriving(Show)

data Cmd = Cmd { cmdName :: FilePath
               , cmdArgs :: [String]
               }

instance Show Cmd where
  show (Cmd name args) = intercalate " " ([name]++args)

makecmd :: ParameterSet -> Cmd
makecmd (ParameterSet _ tm tau bg) = Cmd name args 
  where
    name = exe
    args = [source tm
           ,"1"
           ,show tau
           ,show bg
           ,"74"
           ]

runcmd :: Cmd -> IO ()
runcmd cmd@(Cmd name args) = do
  putStrLn $ "--- Running: "
  putStrLn $ show cmd
  t0<-getCPUTime
  stdout<-readProcess name args [] 
  putStrLn stdout
  t1<-getCPUTime
  let diff = (fromIntegral (t1 - t0)) / (10^12)
  printf "Execution Time: %0.3f sec\n" (diff :: Double)
  putStrLn "---"

select :: Int -> [a] -> [a]
select workid
  | workid>=0 = take batch . drop (batch*workid)
  | otherwise = \x->[]

recommendation :: (RealFrac a) => a->String
recommendation workItemCount = 
  printf "qsub -t 1-%d -pe batch 16 -l short=true -N pstack-params -j y -o /dev/null -b y -cwd -V 'src/Main ${SGE_TASK_ID} > output.${SGE_TASK_ID}'" (1+(ceiling workItemCount)::Int)

main = do
  args <- getArgs
  if (length args==0)
    then 
      --putStrLn $ show $ (fromIntegral (length params) / fromIntegral (batch))
      putStrLn $ recommendation $ (fromIntegral (length params) / fromIntegral (batch)) 
    else do
      mapM (putStrLn . show) $ select (workid args) $ fmap makecmd params'
      mapM runcmd $ select (workid args) $ fmap makecmd params'
      return ()
  
  where  
    indexedparams=zip [0..] params
    params' = [ParameterSet id tm tau bg | (id,(tm,tau,bg))<-indexedparams]
    
    -- input has to use 1-based indexing because of qsub restrictions around array jobs
    -- so: convert here to 0 based indexing
    workid :: [String] -> Int
    workid args = (read $ head args)-1

    
