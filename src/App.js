import { useState } from 'react';
import Navbar from './components/Navbar';
import LandingPage from './pages/LandingPage';
import AnalysisApp from './pages/AnalysisApp';

function App() {
  const [currentPage, setCurrentPage] = useState('landing');
  const goToApp = () => setCurrentPage('app');
  const goToLanding = () => setCurrentPage('landing');

  return (
    <div style={{
      minHeight: '100vh',
      background: '#0a0f1a',
      color: '#e8f0fe',
      fontFamily: "'Segoe UI', sans-serif"
    }}>
      <Navbar onLaunchApp={goToApp} onGoHome={goToLanding} />

      {currentPage === 'landing' && <LandingPage onLaunchApp={goToApp} />}
      {currentPage === 'app' && <AnalysisApp onBack={goToLanding} />}
    </div>
  );
}

export default App;