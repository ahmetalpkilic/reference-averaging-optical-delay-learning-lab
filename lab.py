"""Independent reference averaging; dimensional lag matrices have units s²."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from shared_reference import sequences, lag_batch, score, variance_components, reference_differences
from model import SCALE, threshold, weights, mean_covariance
from whitening import calibration
from optics import image, fit_center, ORIGIN, VELOCITY, READ_SD

SEED, TRAIN, BATCHES, EVENTS, N, M = 20261023, 999, 32, 200, 96, 32
LEVELS, K, TAU = (1, 4, 16), 8, .12


def average_lags(traces):
    """Average separately centered trace lags, NOT lag of the averaged trace.

    Input (independent references, time samples, channels). Each reference
    estimates the same population lag; averaging preserves that expectation.
    """
    x = np.asarray(traces)
    if x.ndim != 3 or x.shape[0] < 1 or x.shape[1] < 2 or x.shape[2] != 2:
        raise ValueError('Expected nonempty (references, time>=2, 2) traces')
    if not np.isfinite(x).all():
        raise ValueError('Reference samples must be finite')
    return lag_batch(x).mean(axis=0)


def nested_references(rng, count):
    """Independent reset traces; nested means pair the comparison across L."""
    total = np.zeros((count, 2, 2)); result = {}
    for ell in range(1, max(LEVELS)+1):
        total += lag_batch(sequences(rng, count, N, 0))
        if ell in LEVELS:
            result[ell] = total.copy()/ell
    return result


def covariance_check():
    """Independent Monte Carlo validation of covariance scaling, in s^4."""
    refs = nested_references(np.random.default_rng(SEED+1), 6000)
    covs = {ell: np.cov(ref[:, [0, 0, 1], [0, 1, 1]], rowvar=False)
            for ell, ref in refs.items()}
    baseline = np.trace(covs[1])
    return [dict(references=ell, expected_ratio=1/ell,
                 measured_trace_ratio=float(np.trace(covs[ell])/baseline),
                 covariance_s4=covs[ell].tolist()) for ell in LEVELS]


def run():
    cov_rng, ref_rng, train_ref_rng, train_pilot_rng, test_rng, clock_rng, camera_rng = [
        np.random.default_rng(s) for s in np.random.SeedSequence(SEED).spawn(7)]
    w0, w1 = weights(M,0), weights(M,60)
    rows=[]
    for batch in range(BATCHES):
        c=calibration(cov_rng,K)  # Same calibration for all L and both policies.
        refs=nested_references(ref_rng,1)
        train_refs=nested_references(train_ref_rng,TRAIN)
        train_pilots=lag_batch(sequences(train_pilot_rng,TRAIN,N,0))
        cuts={}
        for ell in LEVELS:
            # Matched marginal training: reference bank size ALSO equals L.
            # Conditional training freezes the actual deployed reference mean.
            for policy,reference in [('cov-only',train_refs[ell]),('cov+ref',refs[ell])]:
                cuts[ell,policy]=float(threshold(score(train_pilots-reference,c)[:,None])[0])
        for angle in (0,60):
            pilots=lag_batch(sequences(test_rng,EVENTS,N,angle))
            means=sequences(clock_rng,EVENTS,M,angle).mean(axis=1)
            camera=image(ORIGIN+VELOCITY*TAU)+camera_rng.normal(0,READ_SD,(EVENTS,33))
            optical=np.array([(fit_center(y)-ORIGIN)/VELOCITY-TAU for y in camera])
            error=np.column_stack((optical-means@w0,optical-means@w1))*1000
            for ell in LEVELS:
                scores=score(pilots-refs[ell],c)
                for policy in ('cov-only','cov+ref'):
                    alarm=scores>cuts[ell,policy]
                    selected=np.where(alarm[:,None],w1,w0)
                    rows.append(dict(batch=batch,angle_deg=angle,references=ell,policy=policy,
                        covariance_s2=c.tolist(),reference_lag_s2=refs[ell][0].tolist(),
                        threshold=cuts[ell,policy],alarm_rate=float(alarm.mean()),
                        mse_ms2=float(np.mean(np.where(alarm,error[:,1],error[:,0])**2)),
                        clock_mse_ms2=float(1e6*np.einsum('ni,ij,nj->n',selected,mean_covariance(M,angle),selected).mean())))
    summaries=[]
    for angle in (0,60):
        for ell in LEVELS:
            for policy in ('cov-only','cov+ref'):
                sub=[r for r in rows if (r['angle_deg'],r['references'],r['policy'])==(angle,ell,policy)]
                rates=np.array([r['alarm_rate'] for r in sub])
                summaries.append(dict(angle_deg=angle,references=ell,policy=policy,
                    mean_alarm_percent=float(100*rates.mean()),
                    batch_alarm_percent_q10_q90=(100*np.quantile(rates,[.1,.9])).tolist(),
                    alarm_cluster_mcse_percent=float(100*rates.std(ddof=1)/np.sqrt(BATCHES)),
                    pooled_delay_rmse_ms=float(np.sqrt(np.mean([r['mse_ms2'] for r in sub]))),
                    mean_clock_mse_ms2=float(np.mean([r['clock_mse_ms2'] for r in sub])),
                    **variance_components(rates,EVENTS)))
    return dict(seed=SEED,batches=BATCHES,events_per_batch=EVENTS,training_scores=TRAIN,
                calibration_size=K,pilot_length=N,final_repeats=M,summary=summaries,
                reference_covariance_check=covariance_check(),batches_detail=rows)


def plot(result,path):
    fig,axes=plt.subplots(1,3,figsize=(13,4),layout='constrained')
    for ax,angle,key,title in [(axes[0],0,'alarm_rate','Null alarms (%)'),
                              (axes[1],60,'alarm_rate','60° detection (%)'),
                              (axes[2],60,'mse_ms2','60° pooled delay RMSE (ms)')]:
        for j,ell in enumerate(LEVELS):
            for policy,offset,color in [('cov-only',-.12,'#386cb0'),('cov+ref',.12,'#e47920')]:
                v=np.array([r[key] for r in result['batches_detail']
                    if (r['angle_deg'],r['references'],r['policy'])==(angle,ell,policy)])
                center=np.sqrt(v.mean()) if key=='mse_ms2' else 100*v.mean()
                values=np.sqrt(v) if key=='mse_ms2' else 100*v
                ax.scatter(np.full(len(v),j+offset),values,s=10,alpha=.35,color=color,
                           label=policy if j==0 else None)
                ax.plot([j+offset-.09,j+offset+.09],[center]*2,color=color,lw=3)
        ax.set_xticks(range(3),['L=1','L=4','L=16']);ax.set_title(title);ax.grid(axis='y',alpha=.2)
    axes[0].axhline(5,color='gray',ls=':');axes[0].legend(fontsize=8)
    fig.suptitle('Independent reference averaging; paired batches, fixed K=8 and pilot n=96')
    fig.savefig(path,dpi=160);plt.close(fig)

if __name__=='__main__':
    out=Path('examples');out.mkdir(exist_ok=True)
    result=run();(out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    plot(result,out/'averaging.png')
    for r in result['summary']:
        print(f"angle={r['angle_deg']:2} L={r['references']:2} {r['policy']:8} alarm={r['mean_alarm_percent']:5.2f}% RMSE={r['pooled_delay_rmse_ms']:.3f} ms")
    print('Reference covariance:',result['reference_covariance_check'])
